#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "parser.c"

// Polymorphism in C

typedef struct Pixel
{
    double r;    //create a pixel struct like discussed in class
    double g;
    double b;
}Pixel;

typedef struct Image
{   //image struct that keeps track of the data in the input file
    int width;
    int height;
    int color;
    unsigned char *data;
} Image;

Image* ray_casting(char*, int, int, Object**);
int write_image(char*, char*, Image*);
double frad(int, double*, Object**);
double fang(int, double*, Object**);
double* diffuse(int, int, double*, double*, Object**);
double* specular(int, int, double, double*, double*, Object**);
double clamp(double);
double* recursive_shoot(int, double*, double*, Object**, int, int);

int main(int argc, char **argv) 
{
    if(argc > 5)
	{
        fprintf(stderr, "Error: Too many arguments!\n");
        return(1);
    }

    char *width = argv[1];
    char *height = argv[2];
    char *input_filename = argv[3];
    char *output_filename = argv[4];

    Object** objects = malloc(sizeof(Object*) * 128);

    int w = atoi(width);
    int h = atoi(height);
    if(w <= 0)
	{
        fprintf(stderr, "Error: width has to be a positive number!\n");
        return 1;
    }
    if(h <= 0)
	{
        fprintf(stderr, "Error: height has to be a positive number!\n");
        return 1;
    }

    read_scene(input_filename, objects);
    
    Image* buffer = ray_casting(input_filename, w, h, objects);
    
    buffer->width = w;
    buffer->height = h;
    
    write_image("P6", output_filename, buffer);
    
	return 0;
}

static inline double sqr(double v) 
{
  return v * v;
}


static inline void normalize(double* v) 
{
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}


double sphere_intersection(double* Ro, double* Rd, double* Center, double r)
{
   	double a = sqr(Rd[0]) + sqr(Rd[1]) + sqr(Rd[2]);
  	double b = 2 * (Ro[0] * Rd[0] + Ro[1] * Rd[1] + Ro[2] * Rd[2] - Rd[0] * Center[0] - Rd[1] * Center[1] - Rd[2] * Center[2]);
  	double c = sqr(Ro[0]) + sqr(Ro[1]) + sqr(Ro[2]) + sqr(Center[0]) + sqr(Center[1]) + sqr(Center[2]) - 2 * (Ro[0] * Center[0] + Ro[1]*Center[1] + Ro[2]*Center[2]) - sqr(r);
    double d = sqr(b) - 4 * a * c;
    if(d < 0)
	{
		return -1;
    }
    d = sqrt(d);
    
    double t1 = (-b - d) / (2 * a);
    double t2 = (-b + d) / (2 * a);
    if(t1 > 0)
	{
        return t1;
    }
    if(t2 > 0)
	{
        return t2;
    }
    
    return -1;
}

double plane_intersection(double* Ro, double* Rd, double* position, double* normal)
{
    double t = - (normal[0] * Ro[0] + normal[1] * Ro[1] + normal[2] * Ro[2] - normal[0] * position[0] - normal[1] * position[1] - normal[2] * position[2]) / (normal[0] * Rd[0] + normal[1] * Rd[1] + normal[2] * Rd[2]);
    if (t > 0)
	{
        return t;
    }
    
    return -1;
}

int write_data(char ppmnum, FILE *output_file, Image* buffer)
{
    if(ppmnum == '6')
	{
        fwrite(buffer->data, sizeof(Pixel), buffer->width * buffer->height, output_file);
        
        printf("Json file converted to PPM File!\n");
        return(0);
    }
    else if(ppmnum == '3')
	{
        int i, j;
        for(i = 0; i < buffer->height; i++)
		{
            for(j = 0; j < buffer->width; j++)
			{
                fprintf(output_file, "%d %d %d ", buffer->data[i * buffer->width * 3 + j * 3], buffer->data[i * buffer->width * 3 + j * 3 + 1], buffer->data[i * buffer->width * 3 + 2]);
            }
            
            fprintf(output_file, "\n");
        }
        
        printf("The file saved successfully!\n");
        return(0);
        
    }
    else
	{
        fprintf(stderr, "Error: incorrect ppm version\n");
        return(1);
    }
}

int write_image(char* outnumber, char* output_filename, Image* buffer)
{
    int width = buffer->width;
    int height = buffer->height;
    
    char ppmNum = outnumber[1];
    
    FILE *fh = fopen(output_filename, "wb");
    if(fh == NULL)
	{
        fprintf(stderr, "Error: cannot open the file\n");
        return(1);
    }

    char *comment = "# output.ppm";

    fprintf(fh, "P%c\n%s\n%d %d\n%d\n", ppmNum, comment, width, height, 255);
    
    write_data(ppmNum, fh, buffer);
    
    fclose(fh);
    return(0);
}

double* intersect(double* Ro, double* Rd, int object_num, Object** objects)
{
    int closest_object_num = -1;
    
    double closest_t = INFINITY;
    
    int i;
    
    double t;

    for(i = 0; i < object_num; i++)
	{
        if(objects[i]->type == 1)
		{
            t = sphere_intersection(Ro, Rd, objects[i]->sphere.position, objects[i]->sphere.radius);
            if(t)
			{
                if(t > 0 && t <= closest_t)
				{
                    closest_t = t;
                    closest_object_num = i;
                }
            }
            else
			{
                fprintf(stderr, "Error: cannot find distance\n");
                exit(1);
            }
        }
        else if(objects[i]->type == 2)
		{
            t = plane_intersection(Ro, Rd, objects[i]->plane.position, objects[i]->plane.normal);
            if(t)
			{
                if(t > 0 && t <= closest_t)
				{
                    closest_t = t;
                    closest_object_num = i;
                }
            }
            else
			{
                fprintf(stderr, "Error: cannot find distance\n");
                exit (1);
            }
        }
    }
    double* result_v;
    result_v = malloc(sizeof(double) * 2);
    result_v[0] = (double) closest_object_num;
    result_v[1] = closest_t;
    
    return result_v;
}

Image* ray_casting(char* filename, int w, int h, Object** objects)
{
    Image* buffer = (Image*)malloc(sizeof(Image));
    if(objects[0] == NULL)
	{
        fprintf(stderr, "Error: no object found");
        exit(1);
    }

    int camera_found = 0;
    
    double width;
    double height;
    
    int i;
    for(i = 0; objects[i] != 0; i += 1)
	{
        if(objects[i]->type == 0)
		{
            camera_found = 1;
            
            width = objects[i]->camera.width;
            height = objects[i]->camera.height;
            if(width<=0 || height<=0)
			{
                fprintf(stderr, "Error: invalid size for camera/n");
                exit(1);
            }
        }
    }

    if(camera_found == 0)
	{
        fprintf(stderr, "Error: Camera is not found\n");
        exit(1);
    }

    buffer->data = (unsigned char*)malloc(w * h * sizeof(Pixel));
    
    Pixel *pixel = (Pixel*)malloc(sizeof(Pixel));
    
    if(buffer->data == NULL || buffer == NULL)
	{
        fprintf(stderr, "Error: cannot allocate memory\n");
        exit(1);
    }

    double pixwidth = width / w;
    double pixheight = height / h;
    
    int j, k;
    
    double Ro[3] = {0, 0, 0};

    for(k = 0; k < h; k++)
	{
        int count = (h - k - 1) * w * 3;
        
        double vy = -height / 2 + pixheight * (k + 0.5);
        
        for(j = 0; j < w; j++)
		{
            double vx = -width / 2 + pixwidth * (j + 0.5);
            double Rd[3] = {vx, vy, 1};

            normalize(Rd);
            
            int recursive_depth;
            int inside_sphere;
            
            pixel->r = 0;
        	pixel->g = 0;
        	pixel->b = 0;
            
            double* color = recursive_shoot(i, Rd, Ro, objects, recursive_depth, inside_sphere);
            
            pixel->r = color[0];
			pixel->g = color[1];
			pixel->b = color[2];
        	
            buffer->data[count++] = (unsigned char)255 * clamp(pixel->r);
            buffer->data[count++] = (unsigned char)255 * clamp(pixel->g);
            buffer->data[count++] = (unsigned char)255 * clamp(pixel->b);
        }
    }
    return buffer;
}

double frad(int light_index, double* intersection, Object** objects)
{
	double light_position[3];
	double radial[3];
	
	double distance;
	double calc;
	
	light_position[0] = objects[light_index]->light.position[0];
	light_position[1] = objects[light_index]->light.position[1];
	light_position[2] = objects[light_index]->light.position[2];
	
	distance = sqr(light_position[0]-intersection[0]) + sqr(light_position[1]-intersection[1]) + sqr(light_position[2]-intersection[2]);	
	distance = sqrt(distance);
	
	radial[0] = objects[light_index]->light.radial_a0;
	radial[1] = objects[light_index]->light.radial_a1;
	radial[2] = objects[light_index]->light.radial_a2;
	
	if(distance == INFINITY)
	{
		return 1;
	}
	else
	{
		calc = 1 / (radial[2] * sqr(distance) + radial[1] * distance + radial[0]);
		return calc;
	}
}

double fang(int light_index, double* intersection, Object** objects)
{
	double vl[3];
	double vo[3];
	
	double angular;
	
	vl[0] = objects[light_index]->light.direction[0];
	vl[1] = objects[light_index]->light.direction[1];
	vl[2] = objects[light_index]->light.direction[2];
	normalize(vl);
	
	vo[0] = intersection[0] - objects[light_index]->light.position[0];
	vo[1] = intersection[1] - objects[light_index]->light.position[1];
	vo[2] = intersection[2] - objects[light_index]->light.position[2];
	normalize(vo);
	
	double cosa = vl[0] * vo[0] + vl[1] * vo[1] + vl[2] * vo[2];
	if(objects[light_index]->light.angular_a0 != 0)
	{
		angular = objects[light_index]->light.angular_a0;
	}
	else
	{
		return 1;
	}
	
	if(cos(objects[light_index]->light.theta) > cosa)
	{
		return 0;
	}
	else
	{
		return pow(cosa, angular);
	}
}

double* diffuse(int object_index, int light_index, double* N, double* L, Object** objects){
	double val = N[0] * L[0] + N[1] * L[1] + N[2] * L[2];
	
	double* calc;
	double* k;
	
	calc = malloc(sizeof(double) * 3);
	k = malloc(sizeof(double) * 3);
	
	if (val <= 0)
	{
		calc[0] = 0;
		calc[1] = 0;
		calc[2] = 0;
	}
	else 
	{
		if (objects[object_index]->type == 1)
		{
			k[0] = objects[object_index]->sphere.diffuse_color[0] * objects[light_index]->light.color[0];
			k[1] = objects[object_index]->sphere.diffuse_color[1] * objects[light_index]->light.color[1];
			k[2] = objects[object_index]->sphere.diffuse_color[2] * objects[light_index]->light.color[2];
			
			calc[0] = k[0] * val;
			calc[1] = k[1] * val;
			calc[2] = k[2] * val;
		}
		else if (objects[object_index]->type == 2)
		{
			k[0] = objects[object_index]->plane.diffuse_color[0] * objects[light_index]->light.color[0];
			k[1] = objects[object_index]->plane.diffuse_color[1] * objects[light_index]->light.color[1];
			k[2] = objects[object_index]->plane.diffuse_color[2] * objects[light_index]->light.color[2];
			
			calc[0] = k[0] * val;
			calc[1] = k[1] * val;
			calc[2] = k[2] * val;
		}
	}
	return calc;
}

double* specular(int object_index, int light_index, double val, double* V, double* R, Object** objects)
{
	double val2 = V[0] * R[0] + V[1] * R[1] + V[2] * R[2];
	
	double* calc;
	double* k;
	
	calc = malloc(sizeof(double) * 3);
	k = malloc(sizeof(double) * 3);
	
	if (val <= 0 || val2 <= 0)
	{
		calc[0] = 0;
		calc[1] = 0;
		calc[2] = 0;
	}
	else 
	{
		if (objects[object_index]->type == 1)
		{
			k[0] = objects[object_index]->sphere.specular_color[0] * objects[light_index]->light.color[0];
			k[1] = objects[object_index]->sphere.specular_color[1] * objects[light_index]->light.color[1];
			k[2] = objects[object_index]->sphere.specular_color[2] * objects[light_index]->light.color[2];
			
			if(objects[object_index]->light.ns == 0)
			{
				calc[0] = k[0] * pow(val2, 20);
				calc[1] = k[1] * pow(val2, 20);
				calc[2] = k[2] * pow(val2, 20);
			}
			else
			{
				calc[0] = k[0] * pow(val2, objects[object_index]->light.ns);
				calc[1] = k[1] * pow(val2, objects[object_index]->light.ns);
				calc[2] = k[2] * pow(val2, objects[object_index]->light.ns);
			}
		}
		else if (objects[object_index]->type == 2)
		{
			k[0] = objects[object_index]->plane.specular_color[0] * objects[light_index]->light.color[0];
			k[1] = objects[object_index]->plane.specular_color[1] * objects[light_index]->light.color[1];
			k[2] = objects[object_index]->plane.specular_color[2] * objects[light_index]->light.color[2];
		
			if(objects[object_index]->light.ns == 0)
			{
				calc[0] = k[0] * pow(val2, 20);
				calc[1] = k[1] * pow(val2, 20);
				calc[2] = k[2] * pow(val2, 20);
			}
			else
			{
				calc[0] = k[0] * pow(val2, objects[object_index]->light.ns);
				calc[1] = k[1] * pow(val2, objects[object_index]->light.ns);
				calc[2] = k[2] * pow(val2, objects[object_index]->light.ns);
			}
		}
	}
	return calc;
}

double clamp(double number)
{
	if (number <= 0)
	{
		return 0;
	}
	else if (number >= 1)
	{
		return 1;
	}
	else {
		return number;
	}
}

double* recursive_shoot(int object_num, double* Rd, double* Ro, Object** objects, int recursive_depth, int inside_sphere) 
{
	double* color;
	color = malloc(sizeof(double) * 3);
	color[0] = 0;
	color[1] = 0;
	color[2] = 0;
	
	double reflectivity;
	double refractivity;
	double ior;
	
	if (recursive_depth > 7) 
	{
		return color;
	}
	
	double* inter;
	inter = malloc(sizeof(double) * 2);
	inter = intersect(Ro, Rd, object_num, objects);
	
	int intersection = (int)inter[0];
	
	double closest_t = inter[1];
	
	int has_shadow = 0;
	
	if (intersection >= 0) 
	{
		double Ron[3];
		double N[3];
		double V[3];
			
		V[0] = Rd[0];
		V[1] = Rd[1];
		V[2] = Rd[2];
		
		Ron[0] = closest_t * Rd[0] + Ro[0];
		Ron[1] = closest_t * Rd[1] + Ro[1];
		Ron[2] = closest_t * Rd[2] + Ro[2];
		
		if (objects[intersection]->type == 1) 
		{
			N[0] = Ron[0] - objects[intersection]->sphere.position[0];
			N[1] = Ron[1] - objects[intersection]->sphere.position[1];
			N[2] = Ron[2] - objects[intersection]->sphere.position[2]
			;
			reflectivity = objects[intersection]->sphere.reflectivity;
			refractivity = objects[intersection]->sphere.refractivity;
			
			ior = objects[intersection]->sphere.ior;
		}
		else if (objects[intersection]->type == 2) 
		{
			N[0] = objects[intersection]->plane.normal[0];
			N[1] = objects[intersection]->plane.normal[1];
			N[2] = objects[intersection]->plane.normal[2];
			
			reflectivity = objects[intersection]->plane.reflectivity;
			refractivity = objects[intersection]->plane.refractivity;
			
			ior = objects[intersection]->plane.ior;
		}
		
		normalize(N);
		
		double intersect_position[3];
		intersect_position[0] = Ron[0];
		intersect_position[1] = Ron[1];
		intersect_position[2] = Ron[2];
		
		double L[3];
		double R[3];
		double Rdn[3];
		
		int n;
		for (n = 0; objects[n] != 0; n++) 
		{
			if (objects[n]->type == 3) 
			{
				Rdn[0] = objects[n]->light.position[0] - Ron[0];
				Rdn[1] = objects[n]->light.position[1] - Ron[1];
				Rdn[2] = objects[n]->light.position[2] - Ron[2];
				
				double a;
				
				double light_distance = sqrt(sqr(Rdn[0]) + sqr(Rdn[1]) + sqr(Rdn[2]));
				a = light_distance;
				
				normalize(Rdn);
				
				int b;
				for (b = 0; objects[b] != 0; b++) 
				{
					if (b != intersection) 
					{
						if (objects[b]->type == 1) 
						{
							a = sphere_intersection(Ron, Rdn, objects[b]->sphere.position, objects[b]->sphere.radius);
							if (a > 0 && a < light_distance) 
							{
								has_shadow = 1;
							}
						}
						else if (objects[b]->type == 2) 
						{
							a = plane_intersection(Ron, Rdn, objects[b]->plane.position, objects[b]->plane.normal);
							if (a > 0 && a < light_distance) 
							{
								has_shadow = 1;
							}
						}
					}
				}
				if (has_shadow == 0) 
				{
					L[0] = Rdn[0];
					L[1] = Rdn[1];
					L[2] = Rdn[2];
					normalize(L);
				
					double NL = N[0] * L[0] + N[1] * L[1] + N[2] * L[2];
					R[0] = -2 * NL*N[0] + L[0];
					R[1] = -2 * NL*N[1] + L[1];
					R[2] = -2 * NL*N[2] + L[2];
					
					double* diff;
					double* spec;
					diff = malloc(sizeof(double) * 3);
					spec = malloc(sizeof(double) * 3);
					
					double fr, fa;
					fr = frad(n, intersect_position, objects);
					fa = fang(n, intersect_position, objects);
					
					diff = diffuse(intersection, n, N, L, objects);
					spec = specular(intersection, n, NL, V, R, objects);
					
					color[0] += fr*fa*(diff[0] + spec[0]);
					color[1] += fr*fa*(diff[1] + spec[1]);
					color[2] += fr*fa*(diff[2] + spec[2]);

					double new_Ro[3];
					new_Ro[0] = Ron[0];
					new_Ro[1] = Ron[1];
					new_Ro[2] = Ron[2];
					
					double* reflection_color;
					double* refraction_color;
					
					reflection_color = malloc(sizeof(double) * 3);
					refraction_color = malloc(sizeof(double) * 3);
					
					reflection_color[0] = 0;
					reflection_color[1] = 0;
					reflection_color[2] = 0;
					
					refraction_color[0] = 0;
					refraction_color[1] = 0;
					refraction_color[2] = 0;
					
					double new_Rd[3];
					
					if (reflectivity > 0) 
					{
						new_Rd[0] = R[0];
						new_Rd[1] = R[1];
						new_Rd[2] = R[2];
						
						double offset[3] = { 0, 0, 0 };
						offset[0] = new_Rd[0] * 0.0001;
						offset[1] = new_Rd[1] * 0.0001;
						offset[2] = new_Rd[2] * 0.0001;
						
						new_Ro[0] = new_Ro[0] + offset[0];
						new_Ro[1] = new_Ro[1] + offset[1];
						new_Ro[2] = new_Ro[2] + offset[2];
						
						normalize(new_Rd);
						
						reflection_color = recursive_shoot(object_num, new_Rd, new_Ro, objects, recursive_depth + 1, inside_sphere);
					}
					if (refractivity > 0) 
					{
						if (inside_sphere == 1) 
						{
							ior = 1 / ior;
						}
						if (objects[intersection]->type == 1 && inside_sphere == 0) 
						{
							inside_sphere = 1;
						}
						else if (objects[intersection]->type == 1 && inside_sphere == 1) 
						{
							inside_sphere = 0;
						}
						
						double a[3];
						double b[3];
						
						double sinPhi, cosPhi;
				
						a[0] = N[1] * L[2] - N[2] * L[1];
						a[1] = N[2] * L[0] - N[0] * L[2];
						a[2] = N[0] * L[1] - N[1] * L[0];
						normalize(a);
						
						b[0] = a[1] * N[2] - a[2] * N[1];
						b[1] = a[2] * N[0] - a[0] * N[2];
						b[2] = a[0] * N[1] - a[1] * N[0];
						
						sinPhi = ior * (L[0] * b[0] + L[1] * b[1] + L[2] * b[2]);
						cosPhi = sqrt(1 - sqr(sinPhi));
						
						new_Rd[0] = -N[0] * cosPhi + b[0] * sinPhi;
						new_Rd[1] = -N[1] * cosPhi + b[1] * sinPhi;
						new_Rd[2] = -N[2] * cosPhi + b[2] * sinPhi;
						
						double offset[3] = { 0, 0, 0 };
						offset[0] = new_Rd[0] * 0.0001;
						offset[1] = new_Rd[1] * 0.0001;
						offset[2] = new_Rd[2] * 0.0001;
						
						new_Ro[0] = new_Ro[0] + offset[0];
						new_Ro[1] = new_Ro[1] + offset[1];
						new_Ro[2] = new_Ro[2] + offset[2];
						
						normalize(new_Rd);
						
						refraction_color = recursive_shoot(object_num, new_Rd, new_Ro, objects, recursive_depth + 1, inside_sphere);
					}
					color[0] = (1 - reflectivity - refractivity) * color[0] + refraction_color[0] * refractivity + reflection_color[0] * reflectivity;
					color[1] = (1 - reflectivity - refractivity) * color[1] + refraction_color[1] * refractivity + reflection_color[1] * reflectivity;
					color[2] = (1 - reflectivity - refractivity) * color[2] + refraction_color[2] * refractivity + reflection_color[2] * reflectivity;
				}
			}
		}
	}
	return color;
}
