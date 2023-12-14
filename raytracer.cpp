#include <iostream>
#include <stdio.h>
#include <math.h>

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

    vec3(float x = 0.0f, float y = 0.0f, float z = 0.0f) {
        x = x;
        y = y;
        z = z;
    }

    float norm () {
        return sqrtf(powf(x, 2.0f) + powf(y, 2.0) + powf(z, 2.0f));
    }

    vec3 operator*(float scalar) {
        return vec3(x * scalar, y*scalar, z * scalar);
    }

    vec3 operator+(vec3 vector) {
        return vec3(x+vector.x, y+vector.y, z+vector.z);
    }

    vec3 operator-(vec3 vector) {
        return vec3(x-vector.x, y-vector.y, z-vector.z);
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

struct Sphere {
    vec3 scale {1.0, 1.0, 1.0};
    vec3 position;
    float radius;
    float matrix[4][4];

    Sphere(vec3 scale, vec3 position) {
        matrix[0][0] = scale.x;
        matrix[1][1] = scale.y;
        matrix[2][2] = scale.z;
        matrix[0][3] = position.x;
        matrix[1][3] = position.y;
        matrix[2][3] = position.z;
        matrix[3][3] = 1.0;
    }
};

// struct vec4 {
//     float x,y,z,w;

//     vec4(float x = 0.0, float y = 0.0, float z = 0.0, float w = 0.0) {
//         x = 0.0;
//         y = 0.0;
//         z = 0.0;
//         w = 0.0;
//     }
// };

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
        vec3 pointOfIntersection1;
        vec3 pointOfIntersection2;
    };

    IntersectResult testIntersection(Sphere* sphere) {
        /* 

            Solve ray/sphere intersection quadratic:
            a = |dir|^2
            b = 2dir(startingPoint - sphereCentre)
            c = |startingPoint - sphereCentre|^2 - radius^2
        */
       
        // Unit direction vector
        float a = 1;
        float b = (direction * 2.0f).dot(startingPoint - sphere->position);
        float c = powf((startingPoint - sphere->position).norm(), 2.f) - powf(sphere->radius, 2.f);

        float det = -powf(b, 2.0) - 4 * a * c;

        IntersectResult result;
        if (det >= 0.0f) {
            float root1 = (-b + sqrtf(det)) / (2.0 * a);
            float root2 = (-b - sqrtf(det)) / (2.0 * a);

            result.intersect = true;
            result.pointOfIntersection1 = startingPoint + direction * root1;
            result.pointOfIntersection2 = startingPoint + direction * root2;
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


// This main function is meant only to illustrate how to use the save_imageXX functions.
// You should get rid of this code, and just paste the save_imageXX functions into your
// raytrace.cpp code. 
int main() {
    vec3 sphereScale = { 1.0, 1.0, 1.0 };
    vec3 position = {0.0, 0.0, 0.0};
    Sphere test{sphereScale, position};

	int Width = 128;	// Move these to your setup function. The actual resolution will be
	int Height= 128;	// specified in the input file
    char fname3[20] = "sceneP3.ppm"; //This should be set based on the input file
	char fname6[20] = "sceneP6.ppm"; //This should be set based on the input file
	unsigned char *pixels;
	// This will be your image. Note that pixels[0] is the top left of the image and
	// pixels[3*Width*Height-1] is the bottom right of the image.
    pixels = new unsigned char [3*Width*Height]{0};

    RGB col = {1.0, 0.0, 0.0};
    setPixel(pixels, 24, 50, col, 128, 128);

	save_imageP3(Width, Height, fname3, pixels);
	save_imageP6(Width, Height, fname6, pixels);

    vec3 dir = {1, 1, 1};
    Ray testRay{dir, vec3{0, 0, 0}};

    Ray::IntersectResult intersect = testRay.testIntersection(&test);


    float c = 3;
}

// Trace ray returns rgb colour data
vec3 traceRay(Ray start) {
    return vec3{0.0, 0.0, 0.0};
}
