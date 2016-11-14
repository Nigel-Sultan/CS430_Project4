#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "header.h"

void skip_ws(FILE*);
char* next_string(FILE*);
double next_number(FILE*);
double* next_vector(FILE*);
void expect_c(FILE*, int);
int next_c(FILE*);

int line = 1;

void read_scene(char* filename, Object** objects) {
    int c;
    FILE* json = fopen(filename, "r");

    if (json == NULL) {
        fprintf(stderr, "Error: Could not open file \"%s\"\n", filename);
        exit(1);
    }

    skip_ws(json);

    // Find the beginning of the list
    expect_c(json, '[');

    skip_ws(json);

    // Find the objects

    int i = 0;
    while (1) {
        Object* object = malloc(sizeof(Object));
        objects[i] = object;
        c = fgetc(json);
        if (c == ']') {
            fprintf(stderr, "Error: This is the worst scene file EVER.\n");
            fclose(json);
            exit(1);
        }
        if (c == '{') {
            skip_ws(json);

            // Parse the object
            char* key = next_string(json);
            if (strcmp(key, "type") != 0) {
                fprintf(stderr, "Error: Expected \"type\" key on line number %d.\n", line);
                exit(1);
            }

            skip_ws(json);

            expect_c(json, ':');

            skip_ws(json);

            char* value = next_string(json);
            char* temp = value;

            if (strcmp(value, "camera") == 0) {
                objects[i]->type = 0;
            } else if (strcmp(value, "sphere") == 0) {
                objects[i]->type = 1;
            } else if (strcmp(value, "plane") == 0) {
                objects[i]->type = 2;
            } else if(strcmp(value, "light") == 0){
            	objects[i]->type = 3;
			}else {
                fprintf(stderr, "Error: Unknown type, \"%s\", on line number %d.\n", value, line);
                fclose(json);
                exit(1);
            }

            skip_ws(json);

            while (1) {
                c = next_c(json);
                if (c == '}') {
                    // stop parsing this object
                    break;
                } else if (c == ',') {
                    // read another field
                    skip_ws(json);
                    char* key = next_string(json);
                    skip_ws(json);
                    expect_c(json, ':');
                    skip_ws(json);
                    if ((strcmp(key, "width") == 0) ||
                    (strcmp(key, "height") == 0) ||
                    (strcmp(key, "radius") == 0) ||
					(strcmp(key, "radial-a2") == 0) ||
					(strcmp(key, "radial-a1") == 0) ||
					(strcmp(key, "radial-a0") == 0) ||
					(strcmp(key, "angular-a0") == 0) ||
					(strcmp(key, "theta") == 0) ||
					(strcmp(key, "ns") == 0)) {
                        double value = next_number(json);
                        if (strcmp(key, "width") == 0){
                            objects[i]->camera.width = value;
                        }
                        else if(strcmp(key, "height") == 0){
                            objects[i]->camera.height = value;
                        }
                        else if(strcmp(key, "radius") == 0){
                            objects[i]->sphere.radius = value;
                        }
                        else if(strcmp(key, "radial-a2") == 0){
                        	objects[i]->light.radial_a2 = value;
						}
						else if(strcmp(key, "radial-a1") == 0){
							objects[i]->light.radial_a1 = value;
						}
						else if(strcmp(key, "radial-a0") == 0){
							objects[i]->light.radial_a0 = value;
						}
						else if(strcmp(key, "angular-a0") == 0){
							objects[i]->light.angular_a0 = value;
						}
						else if(strcmp(key, "theta") == 0){
							objects[i]->light.theta = value;
						}
						else{
							objects[i]->light.ns = value;
						}
                    }
                    else if ((strcmp(key, "diffuse_color") == 0) ||
                    		(strcmp(key, "specular_color") == 0) ||
							(strcmp(key, "position") == 0) ||
							(strcmp(key, "direction") == 0) ||
                    		(strcmp(key, "normal") == 0) ||
							(strcmp(key, "color") == 0)) {
                        double* value = next_vector(json);
                        if (strcmp(key, "diffuse_color") == 0) {
                        	if (strcmp(temp, "sphere") == 0){
                        		objects[i]->sphere.diffuse_color[0] = value[0];
                            	objects[i]->sphere.diffuse_color[1] = value[1];
                            	objects[i]->sphere.diffuse_color[2] = value[2];
							}
							else if(strcmp(temp, "plane") == 0){
								objects[i]->plane.diffuse_color[0] = value[0];
                            	objects[i]->plane.diffuse_color[1] = value[1];
                            	objects[i]->plane.diffuse_color[2] = value[2];
							}
                        }
                        else if (strcmp(key, "specular_color") == 0){
                        	if (strcmp(temp, "sphere") == 0){
                        		objects[i]->sphere.specular_color[0] = value[0];
                            	objects[i]->sphere.specular_color[1] = value[1];
                            	objects[i]->sphere.specular_color[2] = value[2];
							}
							else if(strcmp(temp, "plane") == 0){
								objects[i]->plane.specular_color[0] = value[0];
                            	objects[i]->plane.specular_color[1] = value[1];
                            	objects[i]->plane.specular_color[2] = value[2];
							}
						}
                        else if (strcmp(key, "position") == 0) {
                            if (strcmp(temp, "sphere") == 0) {
                                objects[i]->sphere.position[0] = value[0];
                                objects[i]->sphere.position[1] = value[1];
                                objects[i]->sphere.position[2] = value[2];
                            }
                            else if (strcmp(temp, "plane") == 0) {
                            	objects[i]->plane.position[0] = value[0];
                            	objects[i]->plane.position[1] = value[1];
                            	objects[i]->plane.position[2] = value[2];
                            }
                            else if (strcmp(temp, "light") == 0){
                            	objects[i]->light.position[0] = value[0];
                            	objects[i]->light.position[1] = value[1];
                            	objects[i]->light.position[2] = value[2];
							}
                            else {
                                fprintf(stderr, "Error: Unknown type!.\n");
                                fclose(json);
                                exit(1);
                            }
                        }
                        else if (strcmp(key, "direction") == 0){
                        	if (strcmp(temp, "light") == 0){
                        		objects[i]->light.direction[0] = value[0];
                            	objects[i]->light.direction[1] = value[1];
                            	objects[i]->light.direction[2] = value[2];
							}
						}
						else if(strcmp(key, "normal") == 0){
							if (strcmp(temp, "plane") == 0){
								objects[i]->plane.normal[0] = value[0];
                            	objects[i]->plane.normal[1] = value[1];
                            	objects[i]->plane.normal[2] = value[2];
							}
						}
						else if(strcmp(key, "color") == 0){
							if(strcmp(temp, "light") == 0){
								objects[i]->light.color[0] = value[0];
                            	objects[i]->light.color[1] = value[1];
                            	objects[i]->light.color[2] = value[2];
							}
						}
                    }
                    else {
                        fprintf(stderr, "Error: Unkonwn property, %s, on line %d.\n", key, line);
                        fclose(json);
                        exit(1);
                    }
                    skip_ws(json);
                } else {
                    fprintf(stderr, "Error: Unexpected value on line %d\n", line);
                    fclose(json);
                    exit(1);
                }
            }
            skip_ws(json);
            c = next_c(json);
            if (c == ',') {
                skip_ws(json);
            } else if (c == ']') {
                fclose(json);
                objects[i + 1] = NULL;
                return;
            } else {
                fprintf(stderr, "Error: Expecting ',' or ']' on line %d.\n", line);
                fclose(json);
                exit(1);
            }
        }
        i = i + 1;
    }
}

// skip_ws() skips white space in the file.
void skip_ws(FILE* json) {
    int c = next_c(json);
    while (isspace(c)) {
        c = next_c(json);
    }
    ungetc(c, json);
}

// expect_c() checks that the next character is d.  If it is not it emits
// an error.
void expect_c(FILE* json, int d) {
    int c = next_c(json);
    if (c == d) return;
    fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, line);
    exit(1);
}

// next_string() gets the next string from the file handle and emits an error
// if a string can not be obtained.
char* next_string(FILE* json) {
    char buffer[129];
    int c = next_c(json);
    if (c != '"') {
        fprintf(stderr, "Error: Expected string on line %d.\n", line);
        exit(1);
    }
    c = next_c(json);
    int i = 0;
    while (c != '"') {
        if (i >= 128) {
            fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
        exit(1);
        }
        if (c == '\\') {
            fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
            exit(1);
        }
        if (c < 32 || c > 126) {
            fprintf(stderr, "Error: Strings may contain only ascii characters on line %d.\n", line);
            exit(1);
        }
        buffer[i] = c;
        i += 1;
        c = next_c(json);
    }
    buffer[i] = 0;
    return strdup(buffer);
}

// next_c() wraps the getc() function and provides error checking and line
// number maintenance
int next_c(FILE* json) {
    int c = fgetc(json);
#ifdef DEBUG
    printf("next_c: '%c'\n", c);
#endif
    if (c == '\n') {
        line += 1;
    }
    if (c == EOF) {
        fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", line);
        exit(1);
    }
    return c;
}

double next_number(FILE* json) {
    double value;
    int f = fscanf(json, "%lf", &value);
    if(f != 1){
        fprintf(stderr, "Error: float is expected\n");
    }
    return value;
}

double* next_vector(FILE* json) {
    double* v = malloc(3*sizeof(double));
    expect_c(json, '[');
    skip_ws(json);
    v[0] = next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[1] = next_number(json);
    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    v[2] = next_number(json);
    skip_ws(json);
    expect_c(json, ']');
    return v;
}
