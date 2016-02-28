#include "evoCA_adder.h"
#include "params.h"

int main(int argc, char* argv[]) {
    if (argc != 9) {
        cout << "Usage: evoCA_adder_test [bit1] [bit2] [bit3] [bit4] [bit5] [bit6] [bit7] [bit8]" << endl;
        cout << "Program will add the first four bits to the second four bits using the cellular automata specified in test_ca.ca and output the result." << endl;
        exit(-1);
    }

    int i, j, cur_step, rule_counter, rule_num;
    int rule_input_vals[RULE_SIZE];
    int *cur_CA, *next_CA;
    int my_num_rules = pow(2, RULE_SIZE);
    c_automata *test_CA = new c_automata();
    FILE *myFile;

    myFile = fopen("test_ca.ca", "r");
    if (myFile == NULL) {
        cout << "CA file test_ca.ca not found." << endl;
        exit(-1);
    }

    fscanf(myFile, "Fitness: %*d\n");
    fscanf(myFile, "Input1 Pos: %d", &test_CA->input_pos[0]);
    for (i = 1; i < BITSTRING_SIZE; i++) {
        fscanf(myFile, ",%d", &test_CA->input_pos[i]);
    }
    fscanf(myFile, "\n");

    fscanf(myFile, "Input2 Pos: %d", &test_CA->input_pos[BITSTRING_SIZE]);
    for (i = (BITSTRING_SIZE + 1); i < (2 * BITSTRING_SIZE); i++) {
        fscanf(myFile, ",%d", &test_CA->input_pos[i]);
    }
    fscanf(myFile, "\n");

    fscanf(myFile, "Output Pos: %d", &test_CA->output_pos[0]);
    for (i = 1; i < (BITSTRING_SIZE + 1); i++) {
        fscanf(myFile, ",%d", &test_CA->output_pos[i]);
    }
    fscanf(myFile, "\n");

    fscanf(myFile, "NumSets: %d\n", &test_CA->num_rulesets);
    test_CA->num_rules = my_num_rules;

    test_CA->rule = new int[my_num_rules * test_CA->num_rulesets];
    fscanf(myFile, "Rule: %d", &test_CA->rule[0]);
    for (i = 1; i < (my_num_rules * test_CA->num_rulesets); i++) {
        fscanf(myFile, ",%d", &test_CA->rule[i]);
    }
    fscanf(myFile, "\n");

    fscanf(myFile, "Rulesets:\n%d", &test_CA->ca_rule_list[0]);
    for (i = 1; i < CA_SIZE; i++) {
        fscanf(myFile, ",%d", &test_CA->ca_rule_list[i]);
    }
    fscanf(myFile, "\n");

    fscanf(myFile, "Init State:\n%d", &test_CA->initial_ca_state[0]);
    for (i = 1; i < CA_SIZE; i++) {
        fscanf(myFile, ",%d", &test_CA->initial_ca_state[i]);
    }
    fscanf(myFile, "\n");
    fclose(myFile);

    rescan_input_pos(test_CA);

    // initialize CA
    cur_CA = new int[CA_SIZE];
    for (i = 0; i < CA_SIZE; i++) {
        if (test_CA->cell_is_input[i]) {
            for (j = 0; j < (2 * BITSTRING_SIZE); j++) {
                if (i == test_CA->input_pos[j]) {
                    break;
                }
            }
            cur_CA[i] = atoi(argv[j + 1]);
        } else {
            cur_CA[i] = test_CA->initial_ca_state[i];
        }
        if (i == 0) {
            cout << endl;
            cout << "CA steps:" << endl;
            cout << cur_CA[0];
        } else {
            cout << "," << cur_CA[i];
        }
    }

    cout << endl;

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
            rule_num = my_num_rules - 1;
            for (j = 1; j <= RULE_SIZE; j++) {
                rule_num -= rule_input_vals[j - 1] * pow(2, (RULE_SIZE - j));
            }
            next_CA[i] = test_CA->rule[rule_num + (test_CA->ca_rule_list[i] * my_num_rules)];

            if (i == 0) {
                cout << next_CA[0];
            } else {
                cout << "," << next_CA[i];
            }

        }
        cout << endl;
        delete [] cur_CA;
        cur_CA = next_CA;
    }

    // output answer
    cout << endl;
    cout << "CA output:" << endl;
    for (i = 0; i < (BITSTRING_SIZE + 1); i++) {
        cout << cur_CA[test_CA->output_pos[i]];
    }
    cout << endl;

    delete test_CA;
    delete [] cur_CA;
}
