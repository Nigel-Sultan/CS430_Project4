#ifndef HEADER_H_
#define HEADER_H_

typedef struct{
  	int type; // 0 = camera, 1 = sphere, 2 = plane, 3 = light
  	union {
    	struct {
      		double width;
      		double height;
    	} camera;
    	struct {
    		double diffuse_color[3];
    		double specular_color[3];
      		double position[3];
      		double radius;
      		double reflectivity;
      		double refractivity;
      		double ior;
    	} sphere;
    	struct {
    		double diffuse_color[3];
    		double specular_color[3];
      		double position[3];
      		double normal[3];
      		double reflectivity;
      		double refractivity;
      		double ior;
    	} plane;
    	struct {
    		double color[3];
      		double position[3];
      		double direction[3];
      		double radial_a2;
      		double radial_a1;
      		double radial_a0;
      		double angular_a0;
      		double theta;
      		double ns;
    	} light;
  	};
} Object;

#endif // HEADER_H_
