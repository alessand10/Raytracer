#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>

// Output in P6 format, a binary file containing:
// P6
// ncolumns nrows
// Max colour value
// colours in binary format thus unreadable
void save_imageP6(int Width, int Height, char* fname,unsigned char* pixels) {
  FILE *fp;
  const int maxVal=255; 
  
  printf("Saving image %s: %d x %d\n", fname,Width,Height);
  fp = fopen(fname,"wb");
  if (!fp) {
        printf("Unable to open file '%s'\n",fname);
        return;
  }
  fprintf(fp, "P6\n");
  fprintf(fp, "%d %d\n", Width, Height);
  fprintf(fp, "%d\n", maxVal);

  for(int j = 0; j < Height; j++) {
		  fwrite(&pixels[j*Width*3], 3,Width,fp);
  }

  fclose(fp);
}

// Output in P3 format, a text file containing:
// P3
// ncolumns nrows
// Max colour value (for us, and usually 255)
// r1 g1 b1 r2 g2 b2 .....
void save_imageP3(int Width, int Height, char* fname,unsigned char* pixels) {
	FILE *fp;
	const int maxVal=255;
	
	printf("Saving image %s: %d x %d\n", fname,Width,Height);
	fp = fopen(fname,"w");
	if (!fp) {
		printf("Unable to open file '%s'\n",fname);
		return;
	}
	fprintf(fp, "P3\n");
	fprintf(fp, "%d %d\n", Width, Height);
	fprintf(fp, "%d\n", maxVal);
	
	int k = 0 ;
	for(int j = 0; j < Height; j++) {
		
		for( int i = 0 ; i < Width; i++)
		{
			fprintf(fp," %d %d %d", pixels[k],pixels[k+1],pixels[k+2]) ;
			k = k + 3 ;
		}
		fprintf(fp,"\n") ;
	}
	fclose(fp);
}

struct RGB {
    float r, g, b;
};

struct vec3 {
    float x,y,z;

    float norm () {
        return sqrtf(powf(x, 2.0f) + powf(y, 2.0) + powf(z, 2.0f));
    }

    vec3 operator*(float scalar) {
        return vec3{x * scalar, y*scalar, z * scalar};
    }

    vec3 operator+(vec3 vector) {
        return vec3{x+vector.x, y+vector.y, z+vector.z};
    }

    vec3 operator-(vec3 vector) {
        return vec3{x-vector.x, y-vector.y, z-vector.z};
    }

    vec3 unit() {
        float len = norm();
        return vec3{x/len, y/len, z/len};
    }


    /**
     * @brief Dot product of this vector and parameter vector
     * 
     * @param vector 
     * @return float 
     */
    float dot(vec3 vector) {
        return x * vector.x + y * vector.y + z * vector.z;
    }
};

struct vec2 {
    float x,y;
};

// In local space all spheres are unit size
struct Sphere {
    vec3 scale {1.0, 1.0, 1.0};
    vec3 position{0.0f, 0.0f, 0.0f};
    float radius = 1.f;
    float matrix[4][4];

    Sphere(vec3 scale, vec3 posIn) {
        matrix[0][0] = scale.x;
        matrix[1][1] = scale.y;
        matrix[2][2] = scale.z;
        matrix[0][3] = posIn.x;
        matrix[1][3] = posIn.y;
        matrix[2][3] = posIn.z;
        matrix[3][3] = 1.0;
        position = posIn;
    }

    /**
     * @brief Given a (assumed) world space intersection
     * with this sphere, returns the normal vector
     * 
     * @param intersectWS 
     * @return vec3 
     */
    vec3 getNormalFromIntersect(vec3 intersectWS) {
        return (intersectWS - position).unit();
    }

};

struct LightSource {
    vec3 location;
    vec3 intensity;
};


struct Ray {
private:
    // Ensure direction is private so it can only be set by the setter (ensuring unit length)
    vec3 direction;
public:
    vec3 startingPoint;

    Ray(vec3 rayDirection, vec3 rayStartingPoint) {
        direction = rayDirection;
        startingPoint = rayStartingPoint;
    }

    void setDirection(vec3 dir) {
        direction = dir.unit();
    }

    vec3 getDirection() {
        return direction;
    }

    struct IntersectResult {
        bool intersect = false;
        // Intersection in world space
        vec3 closestIntersection;
        // The t value r(t) for the nearest intersection. Used for
        // depth comparisons 
        float closestTValue;
        // Intersection normal in world space
        vec3 intersectionNormal;
    };

    IntersectResult testIntersection(Sphere* sphere) {
        /* 

            Solve ray/sphere intersection quadratic:
            a = |dir|^2
            b = 2dir(startingPoint - sphereCentre)
            c = |startingPoint - sphereCentre|^2 - radius^2
        */
       
        // Unit direction vector
        float radius = 1.0f;
        float a = direction.norm();
        float b = ((startingPoint - sphere->position) * 2.f).dot(direction);
        float c = powf((startingPoint - sphere->position).norm(), 2.f) - powf(radius, 2.f);

        float det = powf(b, 2.0) - 4 * a * c;

        IntersectResult result;
        if (det >= 0.0f) {
            float root1 = (-b + sqrtf(det)) / (2.0 * a);
            float root2 = (-b - sqrtf(det)) / (2.0 * a);

            // Prematurely set to true, only set to false if both roots are negative (ray intersects backwards only)
            result.intersect = true;

            // two t's greater than or equal to zero, choose closest 
            if (root2 > 0.0f && root1 > 0.0f) {
                if (root1 < root2) {
                    result.closestIntersection = startingPoint + direction * root1;
                    result.closestTValue = root1;
                }
                else {
                    result.closestIntersection = startingPoint + direction * root2;
                    result.closestTValue = root2;
                }
                
            }
            else if (root2 > 0.0f) {
                result.closestIntersection = startingPoint + direction * root2;
                result.closestTValue = root2;
            }
            else if (root1 > 0.0f) {
                result.closestIntersection = startingPoint + direction * root1;
                result.closestTValue = root1;
            }
            else {
                result.intersect = false;
            }

            // If an intersection occurred, determine the normal of the sphere
            if (result.intersect) {
                result.intersectionNormal = sphere->getNormalFromIntersect(result.closestIntersection);
            }
        }

        return result;

    }
};

/**
 * @brief Sets the pixel colour
 * 
 * @param image 
 * @param x Top left is 0
 * @param y Top left is 0
 * @param colour Values should be normalized
 * @param width 
 * @param height 
 */
void setPixel(unsigned char* image, int x, int y, RGB colour, int width, int height) {
    int entriesPerRow = width * 3;
    int pixelIndex = y * entriesPerRow + x * 3;
    image[pixelIndex] = colour.r * 255.0;
    image[pixelIndex + 1] = colour.g * 255.0;
    image[pixelIndex + 2] = colour.b * 255.0;
}

float near = 1, left = 1, right = -1, bottom = -1, top = 1;
vec2 res{512, 512};

// For this project, the camera is fixed at the origin
// and looking down the -z axis
const vec3 cameraOrigin{0.f, 0.f, 0.f};

/**
 * @brief Creates a ray corresponding to a single pixel in world space
 * starting from the observer and passing through the centre of the
 * desired pixel.
 * 
 * @note pixel [0,0] is the top left [res.x, res.y] is the bottom right.
 * 
 * @param pixel 
 * @return Ray 
 */
Ray pixelToWorldRay(vec2 pixel) {
    // Determine ray in view space before undoing view transform
    // The ray starts at the observer who is at the origin of the world
    vec3 origin = {0.0f, 0.0f, 0.0f};

    // Determine the centre of the pixel in view space:

    // Confusingly, y-axis is the horizontal (pos y left), and x-axis (pos x up) is the vertical
    float pixDY = (left - right) / res.x;
    float pixDX = (top - bottom) / res.y;

    // We subtract half the length/height of the pixel to get to the centre
    float halfPixDY = pixDY / 2.f, halfPixDX = pixDX / 2.f;

    vec3 pixelCentre{top - (pixDX * pixel.y) - halfPixDX, left - (pixDY * pixel.x) - halfPixDY, -near};

    return Ray(pixelCentre - cameraOrigin, cameraOrigin);
}

std::vector<Sphere> spheres = {};
std::vector<LightSource> lightSources = {};

int main() {
    vec3 sphereScale{ 1.0f, 1.0f, 1.0f};
    vec3 position{0.0f, 0.0f, -5.0f};
    Sphere test{sphereScale, position};
    spheres.push(test);

	int Width = res.x;	// Move these to your setup function. The actual resolution will be
	int Height= res.y;	// specified in the input file
    char fname3[20] = "sceneP3.ppm"; //This should be set based on the input file
	char fname6[20] = "sceneP6.ppm"; //This should be set based on the input file
	unsigned char *pixels;
	// This will be your image. Note that pixels[0] is the top left of the image and
	// pixels[3*Width*Height-1] is the bottom right of the image.
    pixels = new unsigned char [3*Width*Height]{0};

    RGB col = {1.0, 0.0, 0.0};

    vec3 dir = {1, 0.0, 0};
    vec3 loc = {-1, 0.75, 0};
    Ray testRay{dir, loc};
    
    Ray::IntersectResult intersect = testRay.testIntersection(&test);
    
    for (float y = 0 ; y < res.x ; y ++) {
        for (float x = 0 ; x < res.y ; x++) {
            Ray worldRay = pixelToWorldRay({x,y});
            if (worldRay.testIntersection(&test).intersect)
                setPixel(pixels, x, y, {0.5, 0.1, 0.1}, res.x, res.y);
        }
    }

    save_imageP3(Width, Height, fname3, pixels);
	save_imageP6(Width, Height, fname6, pixels);

    return 0;
}

/**
 * @brief Fires a ray that tests for intersection against all spheres
 * in the scene.
 * 
 * @param ray 
 * @return Ray::IntersectResult 
 */
Ray::IntersectResult hitTestAllSpheres(Ray ray) {
    Ray::IntersectResult closestResult;
    // Determine ray intersection result
    for (int i = 0 ; i < spheres.size() ; i++) {
        int closestIntersectIndex = -1;

        // Setting this to 0 is okay as long as we factor in if an intersection has been found^
        float closestT = 0.0f;

        Ray::IntersectResult result = ray.testIntersection(spheres[i]);
        if (result.intersect){
            // If there hasn't been an intersection yet (closestIntersectIndex == -1) or we found
            // a closer intersection (result.closestTValue < closestT)
            if (closestIntersectIndex == -1 || result.closestTValue < closestT) {
                closestIntersectIndex = i;
                closestT = result.closestTValue;
                closestResult = result;
            }
        }
    }
}

/**
 * @brief Compute the phong illumination model for single light source
 * 
 * @param V 
 * @param N 
 * @param L 
 * @return vec3 
 */
vec3 computePhongModel(vec3 V, vec3 N, vec3 L) {

}


// Compute the lighting by firing shadow rays to all light sources,
// if there is an unoccluded hit, then use computePhongModel to find their contribution
vec3 computeLighting(vec3 surfacePoint, vec3 surfaceNormal) {

}

/**
 * @brief Recursively trace light rays.
 * 
 * @param ray The ray being traced
 * @param prevBounceSphereIndex The index of the sphere which this ray has reflected from. -1 if initial ray
 * @param remainingBounces The number of bounces remaining
 * @return vec3 
 */
vec3 traceRay(Ray ray, int remainingBounces) {
    // Each bounce triggers a recursive call

    Ray::IntersectResult closestResult = hitTestAllSpheres(ray);

    // We now have a result for closestResult, which is either the nearest hit, or no hit

    // If an intersection occurred, continue recursion, otherwise return no colour contribution

    if (closestResult.intersect){
        if (remainingBounces > 0) {
            // Ray newBounce 
            // vec3 newRayResult = traceRay(newBounce)
        }

    }
    else {
        return {0.0f, 0.0f, 0.0f}
    }
    // 
}

