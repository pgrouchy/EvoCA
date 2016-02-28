#include <sstream>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <climits>
#include <cfloat>
#include <list>

#include "params.h"

using namespace std;

inline double randfloat() {return (rand())/(RAND_MAX+1.0);}

class training_case {

    public:

    int inputs[2 * BITSTRING_SIZE];
    int outputs[BITSTRING_SIZE + 1];

    training_case(int *my_inputs) {
        int i, carry;

        for (i = 0; i < (2 * BITSTRING_SIZE); i++) {
            inputs[i] = my_inputs[i];
        }

        // perform bitwise addition (code from https://stackoverflow.com/questions/13282825/adding-binary-numbers-in-c)
        carry = 0;
        for(i = 0; i < BITSTRING_SIZE; i++) {
            outputs[i] = ((my_inputs[i] ^ my_inputs[BITSTRING_SIZE + i]) ^ carry);
            carry = ((my_inputs[i] & my_inputs[BITSTRING_SIZE + i]) | (my_inputs[i] & carry)) | (my_inputs[BITSTRING_SIZE + i] & carry);
        }
        outputs[BITSTRING_SIZE] = carry;
    };
};

class c_automata {

    public:

    int *rule;
    int ca_rule_list[CA_SIZE];
    int input_pos[2 * BITSTRING_SIZE];
    int output_pos[BITSTRING_SIZE + 1];
    int initial_ca_state[CA_SIZE];
    int fitness;
    int num_rulesets;
    int num_rules;
    bool cell_is_input[CA_SIZE];

    c_automata(){};

    c_automata(int my_num_rules) {
        int i, j;
        double percent_ones;
        bool found_dupe;

        num_rulesets = (int)((double)MAX_NUM_INIT_RULESETS * randfloat()) + 1;
        num_rules = my_num_rules;
        rule = new int[num_rules * num_rulesets];
        for (i = 0; i < num_rulesets; i++) {
            percent_ones = randfloat();
            for (j = 0; j < num_rules; j++) {
                if (randfloat() < percent_ones) {
                    rule[(i * num_rules) + j] = 1;
                } else {
                    rule[(i * num_rules) + j] = 0;
                }
            }
        }

        for (i = 0; i < CA_SIZE; i++) {
            cell_is_input[i] = false;
        }

        i = 0;
        while (i < (2 * BITSTRING_SIZE)) {
            input_pos[i] = (int)((double)CA_SIZE * randfloat());
            found_dupe = false;
            for (j = 0; j < i; j++) {
                if (input_pos[j] == input_pos[i]) {
                    found_dupe = true;
                    break;
                }
            }
            if (!found_dupe) {
                cell_is_input[input_pos[i]] = true;
                i++;
            }
        }
        i = 0;
        while (i < (BITSTRING_SIZE + 1)) {
            output_pos[i] = (int)((double)CA_SIZE * randfloat());
            found_dupe = false;
            for (j = 0; j < i; j++) {
                if (output_pos[j] == output_pos[i]) {
                    found_dupe = true;
                    break;
                }
            }
            if (!found_dupe) {
                i++;
            }
        }
        percent_ones = randfloat();
        for (i = 0; i < CA_SIZE; i++) {
            ca_rule_list[i] = (int)((double)num_rulesets * randfloat());
            if (randfloat() < percent_ones) {
                initial_ca_state[i] = 1;
            } else {
                initial_ca_state[i] = 0;
            }
        }
    };

    c_automata(c_automata *copy_me) {
        int i;
        num_rulesets = copy_me->num_rulesets;
        num_rules = copy_me->num_rules;
        rule = new int[num_rules * num_rulesets];
        for (i = 0; i < (num_rules * num_rulesets); i++) {
            rule[i] = copy_me->rule[i];
        }
        for (i = 0; i < (2 * BITSTRING_SIZE); i++) {
            input_pos[i] = copy_me->input_pos[i];
        }
        for (i = 0; i < (BITSTRING_SIZE + 1); i++) {
            output_pos[i] = copy_me->output_pos[i];
        }
        for (i = 0; i < CA_SIZE; i++) {
            initial_ca_state[i] = copy_me->initial_ca_state[i];
            ca_rule_list[i] = copy_me->ca_rule_list[i];
            cell_is_input[i] = copy_me->cell_is_input[i];
        }
    };

    ~c_automata() {
        delete [] rule;
    };
};

c_automata* tournamentSel(c_automata* *my_pop) {
    int ids[TOURNAMENT_SIZE];
    int fitness[TOURNAMENT_SIZE];
    int best_fit = INT_MAX;
    int i, best_id, counter = 0;
    bool found_dupe;

    while (counter < TOURNAMENT_SIZE) {
        ids[counter] = (int)(POP_SIZE * randfloat());
        found_dupe = false;
        for (i = 0; i < counter; i++) {
            if (ids[i] == ids[counter]) {
                found_dupe = true;
                break;
            }
        }
        if (!found_dupe) {
            fitness[counter] = my_pop[ids[counter]]->fitness;
            counter++;
        }
    }

    // find and return the best individual in the tournament
    for (i = 0; i < TOURNAMENT_SIZE; i++) {
        if (fitness[i] < best_fit) {
            best_fit = fitness[i];
            best_id = ids[i];
        }
    }

    return my_pop[best_id];
}

void rescan_input_pos(c_automata* my_ca) {
    int i, j;

    for (i = 0; i < CA_SIZE; i++) {
        my_ca->cell_is_input[i] = false;
        for (j = 0; j < (2 * BITSTRING_SIZE); j++) {
            if (i == my_ca->input_pos[j]) {
                my_ca->cell_is_input[i] = true;
                break;
            }
        }
    }
}

void save_ca(c_automata* my_ca, ofstream &output_file) {
    int i;

    output_file << "Fitness: " << my_ca->fitness << endl;
    output_file << "Input1 Pos: " << my_ca->input_pos[0];
    for (i = 1; i < BITSTRING_SIZE; i++) {
        output_file << "," << my_ca->input_pos[i];
    }
    output_file << endl;
    output_file << "Input2 Pos: " << my_ca->input_pos[BITSTRING_SIZE];
    for (i = (BITSTRING_SIZE + 1); i < (2 * BITSTRING_SIZE); i++) {
        output_file << "," << my_ca->input_pos[i];
    }
    output_file << endl;
    output_file << "Output Pos: " << my_ca->output_pos[0];
    for (i = 1; i < (BITSTRING_SIZE + 1); i++) {
        output_file << "," << my_ca->output_pos[i];
    }
    output_file << endl;
    output_file << "NumSets: " << my_ca->num_rulesets << endl;
    output_file << "Rule: " << my_ca->rule[0];
    for (i = 1; i < (my_ca->num_rules * my_ca->num_rulesets); i++) {
        output_file << "," << my_ca->rule[i];
    }
    output_file << endl;
    output_file << "Rulesets:" << endl;
    output_file << my_ca->ca_rule_list[0];
    for (i = 1; i < CA_SIZE; i++) {
        output_file << "," << my_ca->ca_rule_list[i];
    }
    output_file << endl;
    output_file << "Init State:" << endl;
    output_file << my_ca->initial_ca_state[0];
    for (i = 1; i < CA_SIZE; i++) {
        output_file << "," << my_ca->initial_ca_state[i];
    }
    output_file << endl;
}

void save_snapshot(c_automata* island_list[][POP_SIZE], string filename, int cur_gen) {
    int isl, pop_id;
    ofstream output_file;

    output_file.open(filename.c_str(), ios::out | ios::trunc);
    output_file << "Gen: " << cur_gen << endl;

    for (isl = 0; isl < NUM_ISLANDS; isl++) {
        output_file << "Isl: " << isl << endl;

        for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {
            save_ca(island_list[isl][pop_id], output_file);
        }
        output_file << endl;
    }

    output_file.close();
}

void load_snapshot(c_automata* island_list[][POP_SIZE], string filename, int &cur_gen) {
    int i, isl, pop_id;
    FILE *input_file;

    input_file = fopen(filename.c_str(), "r");
    fscanf(input_file, "Gen: %d\n", &cur_gen);
    for (isl = 0; isl < NUM_ISLANDS; isl++) {
        fscanf(input_file, "Isl: %*d\n");
        for (pop_id = 0; pop_id < POP_SIZE; pop_id++) {
            island_list[isl][pop_id] = new c_automata();
            island_list[isl][pop_id]->num_rules = pow(2, RULE_SIZE);
            fscanf(input_file, "Fitness: %d\n", &island_list[isl][pop_id]->fitness);
            fscanf(input_file, "Input1 Pos: %d", &island_list[isl][pop_id]->input_pos[0]);
            for (i = 1; i < BITSTRING_SIZE; i++) {
                fscanf(input_file, ",%d", &island_list[isl][pop_id]->input_pos[i]);
            }
            fscanf(input_file, "\n");
            fscanf(input_file, "Input2 Pos: %d", &island_list[isl][pop_id]->input_pos[BITSTRING_SIZE]);
            for (i = (BITSTRING_SIZE + 1); i < (2 * BITSTRING_SIZE); i++) {
                fscanf(input_file, ",%d", &island_list[isl][pop_id]->input_pos[i]);
            }
            fscanf(input_file, "\n");
            fscanf(input_file, "Output Pos: %d", &island_list[isl][pop_id]->output_pos[0]);
            for (i = 1; i < (BITSTRING_SIZE + 1); i++) {
                fscanf(input_file, ",%d", &island_list[isl][pop_id]->output_pos[i]);
            }
            fscanf(input_file, "\n");
            fscanf(input_file, "NumSets: %d\n", &island_list[isl][pop_id]->num_rulesets);
            island_list[isl][pop_id]->rule = new int[island_list[isl][pop_id]->num_rulesets * island_list[isl][pop_id]->num_rules];
            fscanf(input_file, "Rule: %d", &island_list[isl][pop_id]->rule[0]);
            for (i = 1; i < (island_list[isl][pop_id]->num_rulesets * island_list[isl][pop_id]->num_rules); i++) {
                fscanf(input_file, ",%d", &island_list[isl][pop_id]->rule[i]);
            }
            fscanf(input_file, "\n");
            fscanf(input_file, "Rulesets:\n");
            fscanf(input_file, "%d", &island_list[isl][pop_id]->ca_rule_list[0]);
            for (i = 1; i < CA_SIZE; i++) {
                fscanf(input_file, ",%d", &island_list[isl][pop_id]->ca_rule_list[i]);
            }
            fscanf(input_file, "\n");
            fscanf(input_file, "Init State:\n");
            fscanf(input_file, "%d", &island_list[isl][pop_id]->initial_ca_state[0]);
            for (i = 1; i < CA_SIZE; i++) {
                fscanf(input_file, ",%d", &island_list[isl][pop_id]->initial_ca_state[i]);
            }
            fscanf(input_file, "\n");

            rescan_input_pos(island_list[isl][pop_id]);
        }
        fscanf(input_file, "\n");
    }
    fclose(input_file);
}
