#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "../../vendor/cjson/cJSON.h"

#include "P533.h"

typedef struct {
    double tx_lat;
    double tx_lng;
    double rx_lat;
    double rx_lng;
    int year;
    int month;
    int hour;
    int ssn;
    double frequency;
    double txpower;
    double tx_gain;
    double rx_gain;
    const char* name;
} TestScenario;

double deg_to_rad(double degrees) {
    return degrees * M_PI / 180.0;
}

void print_results_json(const char* name, struct PathData* path) {
    printf("{\n");
    printf("  \"name\": \"%s\",\n", name);
    printf("  \"distance\": %.6f,\n", path->distance);
    printf("  \"muf50\": %.6f,\n", path->MUF50);
    printf("  \"bcr\": %.6f,\n", path->BCR);
    printf("  \"snr\": %.6f\n", path->SNR);
    printf("}\n");
}

TestScenario* load_scenarios_from_json(int* num_scenarios) {
    FILE* file = fopen("../test_config.json", "r");
    if (!file) {
        fprintf(stderr, "Error: Cannot open test_config.json\n");
        return NULL;
    }
    
    // Read file content
    fseek(file, 0, SEEK_END);
    long length = ftell(file);
    fseek(file, 0, SEEK_SET);
    
    char* content = malloc(length + 1);
    fread(content, 1, length, file);
    content[length] = '\0';
    fclose(file);
    
    // Parse JSON
    cJSON* json = cJSON_Parse(content);
    if (!json) {
        fprintf(stderr, "Error: Invalid JSON in test_config.json\n");
        free(content);
        return NULL;
    }
    
    cJSON* general = cJSON_GetObjectItem(json, "general");
    cJSON* scenarios_array = cJSON_GetObjectItem(json, "scenarios");
    
    if (!scenarios_array || !cJSON_IsArray(scenarios_array)) {
        fprintf(stderr, "Error: No scenarios array in JSON\n");
        cJSON_Delete(json);
        free(content);
        return NULL;
    }
    
    int count = cJSON_GetArraySize(scenarios_array);
    TestScenario* scenarios = malloc(count * sizeof(TestScenario));
    
    // Get general config values
    double power = cJSON_GetObjectItem(general, "power")->valuedouble;
    double tx_gain = cJSON_GetObjectItem(general, "tx_gain")->valuedouble;
    double rx_gain = cJSON_GetObjectItem(general, "rx_gain")->valuedouble;
    int year = cJSON_GetObjectItem(general, "year")->valueint;
    int month = cJSON_GetObjectItem(general, "month")->valueint;
    
    for (int i = 0; i < count; i++) {
        cJSON* scenario = cJSON_GetArrayItem(scenarios_array, i);
        
        scenarios[i].tx_lat = cJSON_GetObjectItem(scenario, "tx_lat")->valuedouble;
        scenarios[i].tx_lng = cJSON_GetObjectItem(scenario, "tx_lon")->valuedouble;
        scenarios[i].rx_lat = cJSON_GetObjectItem(scenario, "rx_lat")->valuedouble;
        scenarios[i].rx_lng = cJSON_GetObjectItem(scenario, "rx_lon")->valuedouble;
        scenarios[i].frequency = cJSON_GetObjectItem(scenario, "frequency")->valuedouble;
        scenarios[i].hour = cJSON_GetObjectItem(scenario, "hour")->valueint;
        scenarios[i].ssn = cJSON_GetObjectItem(scenario, "ssn")->valueint;
        scenarios[i].year = year;
        scenarios[i].month = month;
        scenarios[i].txpower = power;
        scenarios[i].tx_gain = tx_gain;
        scenarios[i].rx_gain = rx_gain;
        // Allocate and copy name string
        const char* name_str = cJSON_GetObjectItem(scenario, "name")->valuestring;
        char* name_copy = malloc(strlen(name_str) + 1);
        strcpy(name_copy, name_str);
        scenarios[i].name = name_copy;
    }
    
    *num_scenarios = count;
    cJSON_Delete(json);
    free(content);
    return scenarios;
}

int main() {
    int num_scenarios;
    TestScenario* scenarios = load_scenarios_from_json(&num_scenarios);
    
    if (!scenarios) {
        return 1;
    }
    
    printf("[\n");
    
    for (int i = 0; i < num_scenarios; i++) {
        struct PathData path;
        
        // Initialize the path structure to zero
        memset(&path, 0, sizeof(path));
        
        // Set up basic parameters first
        strcpy(path.name, scenarios[i].name);
        strcpy(path.txname, "TX");
        strcpy(path.rxname, "RX");
        
        path.year = scenarios[i].year;
        path.month = scenarios[i].month;
        path.hour = scenarios[i].hour;
        path.SSN = scenarios[i].ssn;
        path.frequency = scenarios[i].frequency;
        path.txpower = scenarios[i].txpower;
        path.BW = 3000.0; // 3 kHz bandwidth
        path.Modulation = 0; // Analog
        path.SorL = 0; // Short path
        path.SNRXXp = 90; // 90% reliability
        path.SNRr = 20.0; // Required SNR
        path.SIRr = 20.0; // Required SIR
        
        // Set locations (convert degrees to radians)
        path.L_tx.lat = deg_to_rad(scenarios[i].tx_lat);
        path.L_tx.lng = deg_to_rad(scenarios[i].tx_lng);
        path.L_rx.lat = deg_to_rad(scenarios[i].rx_lat);
        path.L_rx.lng = deg_to_rad(scenarios[i].rx_lng);
        
        // Allocate memory for the path
        int alloc_result = AllocatePathMemory(&path);
        if (alloc_result != RTN_ALLOCATEP533OK) {
            fprintf(stderr, "Error: Failed to allocate path memory for %s (code: %d)\n", scenarios[i].name, alloc_result);
            continue;
        }
        
        // Set up isotropic antennas
        IsotropicPattern(&path.A_tx, scenarios[i].tx_gain, 1);
        IsotropicPattern(&path.A_rx, scenarios[i].rx_gain, 1);
        
        // Read ionospheric parameters
        char data_path[256] = "../../crates/p533/static_data/";
        int ion_result = ReadIonParametersBin(path.month, path.foF2, path.M3kF2, data_path, 1);
        if (ion_result != RTN_READIONPARAOK) {
            fprintf(stderr, "Error: Failed to read ionospheric parameters for %s (month %d, code: %d)\n", scenarios[i].name, path.month, ion_result);
            FreePathMemory(&path);
            continue;
        }
        
        // Read P1239 data
        char p1239_path[256] = "../../crates/p533/static_data/";
        int p1239_result = ReadP1239(&path, p1239_path);
        if (p1239_result != RTN_READP1239OK) {
            fprintf(stderr, "Error: Failed to read P1239 data for %s (code: %d)\n", scenarios[i].name, p1239_result);
            FreePathMemory(&path);
            continue;
        }
        
        // Run the P533 calculation
        int result = P533(&path);
        if (result == RTN_P533OK) {
            print_results_json(scenarios[i].name, &path);
            if (i < num_scenarios - 1) {
                printf(",\n");
            }
        } else {
            fprintf(stderr, "Error: P533 calculation failed for %s (code: %d)\n", scenarios[i].name, result);
        }
        
        // Free memory
        FreePathMemory(&path);
    }
    
    printf("\n]\n");
    
    // Free scenarios
    for (int i = 0; i < num_scenarios; i++) {
        free((void*)scenarios[i].name);
    }
    free(scenarios);
    
    return 0;
}
