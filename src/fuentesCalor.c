#include"header.h"
#include<math.h>

float calor_3Flat(float x, float y, float t){
	float g1=0, g2=0, g3=0;
	if (0.1 >= x && x<=0.3 && 0.4>=y && y<=0.6) g1 = 0.2;
	if (0.4 >= x && x<=0.6 && 0.4>=y && y<=0.6) g2 = 0.2;
	if (0.7 >= x && x<=0.9 && 0.4>=y && y<=0.6) g3 = 0.2;
	return g1+g2+g3;
}


float calor_3GaussUp(float x, float y, float t){
	float g1=0, g2=0, g3=0;
	if (t<=1200) g1 = 0.3*exp( -((x-0.25)*(x-0.25) + (y-0.75)*(y-0.75)) / (2*(0.1)*(0.1)) );
	if (t<=900) g2 = 0.3*exp( -((x-0.5)*(x-0.25) + (y-0.75)*(y-0.75)) / (2*(0.1)*(0.1)) );
	g3 = 0.3*exp( -((x-0.25)*(x-0.75) + (y-0.75)*(y-0.75)) / (2*(0.1)*(0.1)) );
	return g1+g2+g3;
}


float calor_4GaussSparced(float x, float y, float t){
	float g1=0, g2=0, g3=0, g4=0;
	if (t < 1200) g1 = 0.2*exp( -((x-0.2)*(x-0.2) + (y-0.8)*(y-0.8)) / (2*(0.1)*(0.1)) );
	if (t > 600) g2 = 0.2*exp( -((x-0.2)*(x-0.2) + (y-0.2)*(y-0.2)) / (2*(0.1)*(0.1)) );
	if ( 900 < t && t < 1600) g3 = 0.2*exp( -((x-0.8)*(x-0.8) + (y-0.8)*(y-0.8)) / (2*(0.1)*(0.1)) );
	if (t < 1450)  g4 = 0.2*exp( -((x-0.8)*(x-0.8) + (y-0.2)*(y-0.2)) / (2*(0.1)*(0.1)) );
	return g1+g2+g3+g4;

}
