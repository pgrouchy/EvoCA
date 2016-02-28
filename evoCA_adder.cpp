#include "evoCA_adder.h"
#include "params.h"

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cout << "Usage: evoCA_adder [rand seed]" << endl;
        exit(-1);
    }
    int my_id = atoi(argv[1]);
    //cout << "RunID: " << my_id << endl;
    int i, j, k, isl, isl_2, best_isl, gen, pop_id, cur_train, cur_step, xover_pos, xover_counter, mut_pos, temp_pos;
    int size_training_set = pow(2, (2 * BITSTRING_SIZE));
    int num_rules = pow(2, RULE_SIZE);
    int inputs[2 * BITSTRING_SIZE], new_inputs[2 * BITSTRING_SIZE];
    int *cur_CA, *next_CA, *temp_all_rules;
    int rule_input_vals[RULE_SIZE];
    int rule_counter, rule_num;
    int total_genome_size_bits;
    int xover_ruleset_size, cur_ruleset;
    int temp_num_rulesets;
    int min_fit_isl[NUM_ISLANDS], min_fit_id_isl[NUM_ISLANDS];
    bool found_dupe, did_input_xover, did_xover, did_input_mut, did_mut, ruleset_is_connected;

    int *ruleset_counters, *temp_del_ids;
    int temp_fitness;
    double avg_fit;
    double mut_rate_per_gene, mut_rand;
    double percent_ones;
    int min_fit_id, min_fit_num_rulesets = 1;
    int min_fit = INT_MAX, cur_min_fit, mig_counter;

    ofstream output_file;

    c_automata* island_list[NUM_ISLANDS][POP_SIZE];
    c_automata *offspring, *parent1, *parent2, *next_ca_pop[POP_SIZE];
    training_case *my_training_set[size_training_set];

    srand(my_id * 123456);
    //srand(0);

    // generate full training set
    //cout << "Generating training set (" << size_training_set << " training cases)...";
    for (i = 0; i < (2 * BITSTRING_SIZE); i++) {
        inputs[i] = 0;
    }

    for (i = 0; i < size_training_set; i++) {
        int carry = 0;
        my_training_set[i] = new training_case(inputs);
        if (i != (size_training_set - 1)) {
            new_inputs[0] = ((inputs[0] ^ 1) ^ carry);
            carry = ((inputs[0] & 1) | (inputs[0] & carry)) | (1 & carry);
            for (j = 1; j < (2 * BITSTRING_SIZE); j++) {
                new_inputs[j] = ((inputs[j] ^ 0) ^ carry);
                carry = ((inputs[j] & 0) | (inputs[j] & carry)) | (0 & carry);
            }
            for (j = 0; j < (2 * BITSTRING_SIZE); j++) {
                inputs[j] = new_inputs[j];
            }
        }
    }
    //cout << "done." << endl;

    // initialize population of CAs
    //cout << "Initializing " << POP_SIZE << " CAs...";
    if (!LOAD_POP) {
        gen = 0;
        for (isl = 0; isl < NUM_ISLANDS; isl++) {
            for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {
                island_list[isl][pop_id] = new c_automata(num_rules);
            }
        }
    } else {
        stringstream ss;
        ss << "ca_pop_snap_" << my_id << ".pop";
        load_snapshot(island_list, ss.str(), gen);
    }
    //cout << "done." << endl;

    while (gen <= NUM_GENS) {

        // set all fitnesses to 0
        for (isl = 0; isl < NUM_ISLANDS; isl++) {
            for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {
                island_list[isl][pop_id]->fitness = 0;
            }
        }

        // take a snapshot of the population
        if ((gen > 0) && ((gen % SAVE_INTERVAL) == 0)) {
            stringstream ss;
            ss << "ca_pop_snap_" << my_id << ".pop";
            save_snapshot(island_list, ss.str(), gen);
        }

        // evaluate population
        for (isl = 0; isl < NUM_ISLANDS; isl++) {
            for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {

                // test over the entire training set
                for (cur_train = 0; cur_train < size_training_set; cur_train++) {
                    // initialize CA
                    cur_CA = new int[CA_SIZE];
                    for (i = 0; i < CA_SIZE; i++) {
                        if (island_list[isl][pop_id]->cell_is_input[i]) {
                            for (j = 0; j < (2 * BITSTRING_SIZE); j++) {
                                if (i == island_list[isl][pop_id]->input_pos[j]) {
                                    break;
                                }
                            }
                            if (j < BITSTRING_SIZE) {
                                cur_CA[i] = my_training_set[cur_train]->inputs[BITSTRING_SIZE - 1 - j];
                            } else {
                                cur_CA[i] = my_training_set[cur_train]->inputs[(3 * BITSTRING_SIZE) - 1 - j];
                            }
                        } else {
                            cur_CA[i] = island_list[isl][pop_id]->initial_ca_state[i];
                        }
                    }

                    // run CA for NUM_CA_STEPS
                    for (cur_step = 0; cur_step < NUM_CA_STEPS; cur_step++) {
                        next_CA = new int[CA_SIZE];
                        for (i = 0; i < CA_SIZE; i++) {
                            rule_counter = 0;
                            for (j = (i - ((RULE_SIZE - 1) / 2)); j <= (i + ((RULE_SIZE - 1) / 2)); j++) {
                                if (j < 0) {
                                    rule_input_vals[rule_counter] = cur_CA[j + CA_SIZE];
                                } else if (j >= CA_SIZE) {
                                    rule_input_vals[rule_counter] = cur_CA[j - CA_SIZE];
                                } else {
                                    rule_input_vals[rule_counter] = cur_CA[j];
                                }
                                rule_counter++;
                            }
                            rule_num = num_rules - 1;
                            for (j = 1; j <= RULE_SIZE; j++) {
                                rule_num -= rule_input_vals[j - 1] * pow(2, (RULE_SIZE - j));
                            }
                            next_CA[i] = island_list[isl][pop_id]->rule[rule_num + (island_list[isl][pop_id]->ca_rule_list[i] * num_rules)];
                        }
                        delete [] cur_CA;
                        cur_CA = next_CA;
                    }

                    // calculate fitness
                    temp_fitness = 0;
                    for (i = 0; i < (BITSTRING_SIZE + 1); i++) {
                        temp_fitness += (cur_CA[island_list[isl][pop_id]->output_pos[i]] ^ my_training_set[cur_train]->outputs[BITSTRING_SIZE - i]);
                    }
                    island_list[isl][pop_id]->fitness += temp_fitness;

                    delete [] cur_CA;
                }
            }
        }

        // output stats
        avg_fit = 0.0;
        cur_min_fit = INT_MAX;
        for (isl = 0; isl < NUM_ISLANDS; isl++) {
            min_fit_isl[isl] = INT_MAX;
            for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {
                //avg_fit += ((double)island_list[isl][pop_id]->fitness / ((double)POP_SIZE * (double)NUM_ISLANDS));
                if (island_list[isl][pop_id]->fitness < cur_min_fit) {
                    cur_min_fit = island_list[isl][pop_id]->fitness;
                    min_fit_id = pop_id;
                    best_isl = isl;
                }
                if (island_list[isl][pop_id]->fitness < min_fit_isl[isl]) {
                    min_fit_isl[isl] = island_list[isl][pop_id]->fitness;
                    min_fit_id_isl[isl] = pop_id;
                }
                if ((island_list[isl][pop_id]->fitness < min_fit) || ((island_list[isl][pop_id]->fitness == min_fit) && (island_list[isl][pop_id]->num_rulesets < min_fit_num_rulesets))) {
                    min_fit = island_list[isl][pop_id]->fitness;
                    min_fit_num_rulesets = island_list[isl][pop_id]->num_rulesets;
                    stringstream ss;
                    ss << "best_CA_" << my_id << ".ca";
                    output_file.open(ss.str().c_str(), ios::out | ios::trunc);
                    save_ca(island_list[isl][pop_id], output_file);
                    output_file.close();
                }
            }
        }

        for (isl = 0; isl < NUM_ISLANDS; isl++) { 
            avg_fit += (min_fit_isl[isl] / (double)NUM_ISLANDS);
        }


        cout << "Gen: " << gen << " avgFit: " << avg_fit << " curBestFit: " << cur_min_fit << " overallBestFit: " << min_fit << " bestIsl: " << best_isl << endl;

        // generate next generation
        if (gen != NUM_GENS) {

            for (isl = 0; isl < NUM_ISLANDS; isl++) {

                // elitism
                next_ca_pop[0] = new c_automata(island_list[isl][min_fit_id_isl[isl]]);

                // new random genomes
                for (pop_id = 1; pop_id < (NUM_RAND_INIT + 1); pop_id++)  {
                    next_ca_pop[pop_id] = new c_automata(num_rules);
                }

                for (pop_id = (NUM_RAND_INIT + 1); pop_id < POP_SIZE; pop_id++) {
                //for (pop_id = 0; pop_id < POP_SIZE; pop_id++) 
                    // selection
                    parent1 = tournamentSel(island_list[isl]);
                    offspring = new c_automata(parent1);

                    did_input_mut = false;
                    did_input_xover = false;
                    did_xover = false;
                    did_mut = false;

                    total_genome_size_bits = (2 * BITSTRING_SIZE) + (BITSTRING_SIZE + 1) + (num_rules * offspring->num_rulesets) + CA_SIZE;
                    if (offspring->num_rulesets > 1) {
                        total_genome_size_bits += CA_SIZE;
                    }

                    // sexual recombination
                    if (randfloat() < (double)XOVER_RATE) {
                        did_xover = true;
                        do {
                            parent2 = tournamentSel(island_list[isl]);
                        } while (parent2 == parent1);

                        xover_pos = (int)((double)(total_genome_size_bits - 1) * randfloat());
                        xover_counter = 0;
                        if (parent1->num_rulesets > parent2->num_rulesets) {
                            xover_ruleset_size = num_rules * parent2->num_rulesets;
                        } else {
                            xover_ruleset_size = num_rules * parent1->num_rulesets;
                        }

                        for (i = 0; i < xover_ruleset_size; i++) { // always use smaller of two parent rulesets for xover
                            if (xover_counter > xover_pos) {
                                offspring->rule[i] = parent2->rule[i];
                            }
                            xover_counter++;
                        }

                        for (i = 0; i < (2 * BITSTRING_SIZE); i++) {
                            if (xover_counter > xover_pos) {
                                found_dupe = false;
                                for (j = i; j >= 0; j--) {
                                    if (offspring->input_pos[j] == parent2->input_pos[i]) {
                                        found_dupe = true;
                                        break;
                                    }
                                }
                                if (!found_dupe) {
                                    offspring->input_pos[i] = parent2->input_pos[i];
                                    did_input_xover = true;
                                }
                            }
                            xover_counter++;
                        }

                        for (i = 0; i < (BITSTRING_SIZE + 1); i++) {
                            if (xover_counter > xover_pos) {
                                found_dupe = false;
                                for (j = i; j >= 0; j--) {
                                    if (offspring->output_pos[j] == parent2->output_pos[i]) {
                                        found_dupe = true;
                                        break;
                                    }
                                }
                                if (!found_dupe) {
                                    offspring->output_pos[i] = parent2->output_pos[i];
                                }
                            }
                            xover_counter++;
                        }

                        if (parent1->num_rulesets > 1) {
                            for (i = 0; i < CA_SIZE; i++) {
                                if (xover_counter > xover_pos) {
                                    if (parent2->ca_rule_list[i] < offspring->num_rulesets) {
                                        offspring->ca_rule_list[i] = parent2->ca_rule_list[i];
                                    }
                                }
                                xover_counter++;
                            }
                        }

                        for (i = 0; i < CA_SIZE; i++) {
                            if (xover_counter > xover_pos) {
                                offspring->initial_ca_state[i] = parent2->initial_ca_state[i];
                            }
                            xover_counter++;
                        }

                    }

                    if (did_input_xover) {
                        // find which cells are inputs
                        rescan_input_pos(offspring);
                    }

                    do {
                        // mutation
                        if (offspring->num_rulesets < MAX_NUM_RULESETS) {
                            mut_rate_per_gene = (double)MUT_RATE / (double)(total_genome_size_bits + 1);
                        } else {
                            mut_rate_per_gene = (double)MUT_RATE / (double)(total_genome_size_bits);
                        }
                        for (i = 0; i < offspring->num_rulesets; i++) {
                            for (j = 0; j < num_rules; j++) {
                                if (randfloat() < mut_rate_per_gene) {
                                    offspring->rule[(i * num_rules) + j] = offspring->rule[(i * num_rules) + j] ^ 1;
                                    if (offspring->num_rulesets == 1) {
                                        did_mut = true;
                                    } else {
                                        ruleset_is_connected = false;
                                        for (k = 0; k < CA_SIZE; k++) {
                                            if (offspring->ca_rule_list[k] == i) {
                                                ruleset_is_connected = true;
                                                break;
                                            }
                                        }
                                        if (ruleset_is_connected) {
                                            did_mut = true;
                                        }
                                    }
                                }
                            }
                        }

                        // add new ruleset mutation
                        if ((offspring->num_rulesets < MAX_NUM_RULESETS) && (randfloat() < mut_rate_per_gene)) {
                            offspring->num_rulesets++;
                            temp_all_rules = new int[num_rules * offspring->num_rulesets];
                            for (i = 0; i < (num_rules * (offspring->num_rulesets - 1)); i++) {
                                temp_all_rules[i] = offspring->rule[i];
                            }
                            percent_ones = randfloat();
                            for (i = 0; i < num_rules; i++) {
                                if (randfloat() < percent_ones) {
                                    temp_all_rules[i + (num_rules * (offspring->num_rulesets - 1))] = 1;
                                } else {
                                    temp_all_rules[i + (num_rules * (offspring->num_rulesets - 1))] = 0;
                                }
                            }
                            delete [] offspring->rule;
                            offspring->rule = temp_all_rules;
                        }

                        // input mutation
                        for (i = 0; i < (2 * BITSTRING_SIZE); i++) {
                            if (randfloat() < mut_rate_per_gene) {
                                mut_rand = randfloat();
                                did_mut = true;
                                did_input_mut = true;
                                if ((CA_SIZE == (2 * BITSTRING_SIZE)) || (mut_rand < (1.0 / 3.0))) {
                                    do {
                                        mut_pos = (int)((double)(2 * BITSTRING_SIZE) * randfloat());
                                    } while (mut_pos == i);
                                    temp_pos = offspring->input_pos[mut_pos];
                                    offspring->input_pos[mut_pos] = offspring->input_pos[i];
                                    offspring->input_pos[i] = temp_pos;
                                } else if (mut_rand < (2.0 / 3.0)) {
                                    temp_pos = offspring->input_pos[i];
                                    do {
                                        offspring->input_pos[i] = (int)((double)CA_SIZE * randfloat());
                                        found_dupe = false;
                                        if (offspring->input_pos[i] == temp_pos) {
                                            found_dupe = true;
                                        } else {
                                            for (j = 0; j < (2 * BITSTRING_SIZE); j++) {
                                                if ((j != i) && (offspring->input_pos[i] == offspring->input_pos[j])) {
                                                    found_dupe = true;
                                                    break;
                                                }
                                            }
                                        }
                                    } while (found_dupe);
                                } else {
                                    if (randfloat() < 0.5) {
                                        do {
                                            offspring->input_pos[i] = offspring->input_pos[i] - 1;
                                            if (offspring->input_pos[i] < 0) {
                                                offspring->input_pos[i] += CA_SIZE;
                                            }
                                            found_dupe = false;
                                            for (j = 0; j < (2 * BITSTRING_SIZE); j++) {
                                                if ((j != i) && (offspring->input_pos[i] == offspring->input_pos[j])) {
                                                    found_dupe = true;
                                                    break;
                                                }
                                            }
                                        } while (found_dupe);
                                    } else {
                                        do {
                                            offspring->input_pos[i] = offspring->input_pos[i] + 1;
                                            if (offspring->input_pos[i] >= CA_SIZE) {
                                                offspring->input_pos[i] -= CA_SIZE;
                                            }
                                            found_dupe = false;
                                            for (j = 0; j < (2 * BITSTRING_SIZE); j++) {
                                                if ((j != i) && (offspring->input_pos[i] == offspring->input_pos[j])) {
                                                    found_dupe = true;
                                                    break;
                                                }
                                            }
                                        } while (found_dupe);
                                    }
                                }
                            }
                        }

                        if (did_input_mut) {
                            // find which cells are inputs
                            rescan_input_pos(offspring);
                        }

                        // output mutation
                        for (i = 0; i < (BITSTRING_SIZE + 1); i++) {
                            if (randfloat() < mut_rate_per_gene) {
                                did_mut = true;
                                mut_rand = randfloat();
                                if (mut_rand < (1.0 / 3.0)) {
                                    temp_pos = offspring->output_pos[i];
                                    do {
                                        offspring->output_pos[i] = (int)((double)CA_SIZE * randfloat());
                                        found_dupe = false;
                                        if (offspring->output_pos[i] == temp_pos) {
                                            found_dupe = true;
                                        } else {
                                            for (j = 0; j < (BITSTRING_SIZE + 1); j++) {
                                                if ((j != i) && (offspring->output_pos[i] == offspring->output_pos[j])) {
                                                    found_dupe = true;
                                                    break;
                                                }
                                            }
                                        }
                                    } while (found_dupe);
                                } else if (mut_rand < (2.0 / 3.0)) {
                                    do {
                                        mut_pos = (int)((double)(BITSTRING_SIZE + 1) * randfloat());
                                    } while (mut_pos == i);
                                    temp_pos = offspring->output_pos[mut_pos];
                                    offspring->output_pos[mut_pos] = offspring->output_pos[i];
                                    offspring->output_pos[i] = temp_pos;
                                } else {
                                    if (randfloat() < 0.5) {
                                        do {
                                            offspring->output_pos[i] = offspring->output_pos[i] - 1;
                                            if (offspring->output_pos[i] < 0) {
                                                offspring->output_pos[i] += CA_SIZE;
                                            }
                                            found_dupe = false;
                                            for (j = 0; j < (BITSTRING_SIZE + 1); j++) {
                                                if ((j != i) && (offspring->output_pos[i] == offspring->output_pos[j])) {
                                                    found_dupe = true;
                                                    break;
                                                }
                                            }
                                        } while (found_dupe);
                                    } else {
                                        do {
                                            offspring->output_pos[i] = offspring->output_pos[i] + 1;
                                            if (offspring->output_pos[i] >= CA_SIZE) {
                                                offspring->output_pos[i] -= CA_SIZE;
                                            }
                                            found_dupe = false;
                                            for (j = 0; j < (BITSTRING_SIZE + 1); j++) {
                                                if ((j != i) && (offspring->output_pos[i] == offspring->output_pos[j])) {
                                                    found_dupe = true;
                                                    break;
                                                }
                                            }
                                        } while (found_dupe);
                                    }
                                }
                            }
                        }

                        // initial grid mutation
                        for (i = 0; i < CA_SIZE; i++) {
                            if (randfloat() < mut_rate_per_gene) {
                                offspring->initial_ca_state[i] = offspring->initial_ca_state[i] ^ 1;
                                if (!offspring->cell_is_input[i]) {
                                    did_mut = true;
                                }
                            }
                        }

                        // rule list mutation
                        if (offspring->num_rulesets > 1) {
                            for (i = 0; i < CA_SIZE; i++) {
                                if (randfloat() < mut_rate_per_gene) {
                                    did_mut = true;
                                    cur_ruleset = offspring->ca_rule_list[i];
                                    do {
                                        offspring->ca_rule_list[i] = (int)((double)offspring->num_rulesets * randfloat());
                                    } while (offspring->ca_rule_list[i] == cur_ruleset);
                                }
                            }
                        }
                    } while ((!did_xover) && (!did_mut));

                    next_ca_pop[pop_id] = offspring;
                }
            
                // remove unused rulesets
                for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {
                    if ((next_ca_pop[pop_id]->num_rulesets > 1) && (randfloat() < RULESET_CLEANUP_RATE)) {

                        temp_num_rulesets = next_ca_pop[pop_id]->num_rulesets;

                        ruleset_counters = new int[temp_num_rulesets];
                        for (i = 0; i < temp_num_rulesets; i++) {
                            ruleset_counters[i] = 0;
                        }
                        for (i = 0; i < CA_SIZE; i++) {
                            ruleset_counters[next_ca_pop[pop_id]->ca_rule_list[i]]++;
                        }

                        rule_counter = 0;
                        for (i = 0; i < temp_num_rulesets; i++) {
                            if (ruleset_counters[i] != 0) {
                                rule_counter++;
                            }
                        }

                        temp_all_rules = new int[num_rules * rule_counter];

                        rule_counter = 0;
                        for (i = 0; i < temp_num_rulesets; i++) {
                            if (ruleset_counters[i] != 0) { // ruleset is still in use, keep it
                                for (j = 0; j < num_rules; j++) {
                                    temp_all_rules[rule_counter] = next_ca_pop[pop_id]->rule[(i * num_rules) + j];
                                    rule_counter++;
                                }
                            } else { // ruleset not copied (deleted)
                                next_ca_pop[pop_id]->num_rulesets--;
                            }
                        }

                        delete [] next_ca_pop[pop_id]->rule;
                        next_ca_pop[pop_id]->rule = temp_all_rules;
                        // update references to remaining rulesets
                        rule_counter = temp_num_rulesets - next_ca_pop[pop_id]->num_rulesets;
                        if (rule_counter > 0) {
                            temp_del_ids = new int[rule_counter];
                            j = 0;
                            for (i = 0; i < temp_num_rulesets; i++) {
                                if (ruleset_counters[i] == 0) {
                                    temp_del_ids[j] = i;
                                    j++;
                                }
                            }
                            for (i = 0; i < CA_SIZE; i++) {
                                cur_ruleset = next_ca_pop[pop_id]->ca_rule_list[i];
                                for (j = 0; j < rule_counter; j++) {
                                    if (cur_ruleset > temp_del_ids[j]) {
                                        next_ca_pop[pop_id]->ca_rule_list[i]--;
                                    }
                                }
                            }
                            delete [] temp_del_ids;
                        }

                        delete [] ruleset_counters;
                    }

                    delete island_list[isl][pop_id];
                    island_list[isl][pop_id] = next_ca_pop[pop_id];
                }
            }

            // migration
            if ((gen > 0) && ((gen % MIG_INTERVAL) == 0)) {
                for (isl = 0; isl < NUM_ISLANDS; isl++) {
                    mig_counter = 0;
                    isl_2 = isl + 1;
                    if (isl_2 == NUM_ISLANDS) {
                        isl_2 = 0;
                    }
                    while (mig_counter < MIG_NUM) {
                        pop_id = (int)((double)POP_SIZE * randfloat());
                        i = (int)((double)POP_SIZE * randfloat());
                        delete island_list[isl][pop_id];
                        island_list[isl][pop_id] = new c_automata(island_list[isl_2][i]);
                        mig_counter++;
                    }
                }
            }
        }
        gen++;
    }

    for (isl = 0; isl < NUM_ISLANDS; isl++) {
        for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {
            delete island_list[isl][pop_id];
        }
    }
    for (i = 0; i < size_training_set; i++) {
        delete my_training_set[i];
    }
}
