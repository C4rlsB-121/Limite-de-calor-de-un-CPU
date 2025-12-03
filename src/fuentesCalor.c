#include"header.h"
#include<math.h>

float calor_1Gauss(float x, float y, float t){
	float valor;
	if (t <= 30 && t < 1200){
		valor = exp(- (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) / (2*0.01) );
	}
	else{
		valor = 0;
	}
	return valor;
}

float calor_1Flat(float x, float y, float t){
	if ((x >= 0.4 && x <= 0.6) && (y >= 0.4 && y <= 0.6)){
		return 1;
	}
	return 0;
}


float calor_3Gauss(float x, float y, float t){
	float v1, v2, v3;
	float v = 0;
	if ( t>120 && t< 1600){
		v1 = valor = exp(- (x-0.25)*(x-0.25) + (y-0.75)*(y-0.75) ) / (2*0.01) );
		v2 = valor = exp(- (x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) / (2*0.01) );
		v3 = valor = exp(- (x-0.75)*(x-0.75) + (y-0.75)*(y-0.75) ) / (2*0.01) );
		v = v1 + v2 + v3;
	}
	return v;
}


float calor_4GaussSparced(float x, float y, float t){
	float v1 = 0, v2 = 0, v3 = 0, v4 = 0;
	float v;
	if (100 < t && t < 1100) v1 = exp(- 0.4*(x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) / (2*0.01) );
	if (400 < t && t < 1000) v2 = exp(- 0.3*(x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) / (2*0.01) );
	if (900 < t && t < 1600) v3 = exp(- 0.6*(x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) / (2*0.01) );
	if (t < 200 && t < 900)  v4 = exp(- 0.7*(x-0.5)*(x-0.5) + (y-0.5)*(y-0.5) ) / (2*0.01) );
	v = v1 + v2 + v3 + v4;
	return v;

}
