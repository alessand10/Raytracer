#include <iostream>
#include <stdio.h>

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

    vec3() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
    }
};

struct Sphere {
    vec3 scale;
    vec3 position;
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

struct vec4 {
    float x,y,z,w;

    vec4() {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 0.0;
    }
};

struct Ray {
    vec4 direction;
    vec3 startingPoint;
    float testIntersection(Sphere* sphere) {
        /* 
            Intersection of a plane sphere is simple, we are
            simply looking for where the distance from the ray to
            the center of the circle is exactly r

            We are solving for a t value such that dist(ray(t), sphereCentre) = r
            ray(t)_x = dir_x * t + start_x
            ray(t)_y = dir_y * t + start_y
            ray(t)_z = dir_z * t + start_z

            length(ray(t) - sphereCentre) = r
            (ray(t)_x - centre_x)^2 + (ray(t)_y - centre_y)^2 + (ray(t)_z - centre_z)^2 = r^2
            ray(t)_x



        */

       
       



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
    Sphere test = Sphere{vec3{1.0f, 1.0f, 1.0f}, vec3{}}

	int Width = 128;	// Move these to your setup function. The actual resolution will be
	int Height= 128;	// specified in the input file
    char fname3[20] = "sceneP3.ppm"; //This should be set based on the input file
	char fname6[20] = "sceneP6.ppm"; //This should be set based on the input file
	unsigned char *pixels;
	// This will be your image. Note that pixels[0] is the top left of the image and
	// pixels[3*Width*Height-1] is the bottom right of the image.
    pixels = new unsigned char [3*Width*Height];

    RGB col = {0.0, 1.0, 0.0};
    setPixel(pixels, 24, 50, col, 128, 128);

	save_imageP3(Width, Height, fname3, pixels);
	save_imageP6(Width, Height, fname6, pixels);
}

// Trace ray returns rgb colour data
vec3 traceRay(Ray start) {

}
