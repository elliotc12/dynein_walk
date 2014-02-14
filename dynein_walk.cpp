#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

const double tail_stalk_length = 10;
const double md_stalk_length = 10;

double randAngle(double range);
double dist(double d, double h, double i, double j);

class Dynein {
	public:
		void set_blx(double d);
		void set_bly(double d);
		void set_bla(double d);
		void set_mla(double d);
		void set_ta(double d);
		void set_mra(double d);
		
		double get_bly();
		double get_blx();
		double get_bla();
		double get_mla();
		double get_ta();
		double get_mra();
		
		void log();
		
	private:
		double tx, ty; //Tail domain
		double mlx, mly; //Motor domain
		double mrx, mry;
		double blx, bly; //Binding domain
		double brx, bry;
};

/*** Accessor Methods ***/

void Dynein::set_blx(double d) {
	blx = d;
}

void Dynein::set_bly(double d) {
	bly = d;
}

void Dynein::set_bla(double d) {
	mlx = blx + md_stalk_length*cos(d);
	mly = bly + md_stalk_length*sin(d);
	set_mla(get_mla());
}

void Dynein::set_mla(double d) {
	tx = mlx + tail_stalk_length*cos(d + get_bla());
	ty = mly + tail_stalk_length*sin(d + get_bla());
	set_ta(get_ta());
}

void Dynein::set_ta(double d) {
	mrx = tx + tail_stalk_length*cos(d + get_mla() + get_bla());
	mry = ty + tail_stalk_length*sin(d + get_mla() + get_bla());
	set_mra(get_mra());
}

void Dynein::set_mra(double d) {
	brx = mrx + md_stalk_length*cos(d + get_ta() + get_mla() + get_bla());
	bry = mry + md_stalk_length*sin(d + get_ta() + get_mla() + get_bla());
}

double Dynein::get_bly() {
	return bly;
}
double Dynein::get_blx(){
	return blx;
}
double Dynein::get_bla() {
	return (mly > bly) ? acos((mlx-blx)/dist(blx,bly,mlx,mly)) : -acos((mlx-blx)/dist(blx,bly,mlx,mly));
}
double Dynein::get_mla() {
	return (ty > mly) ? acos((tx-mlx)/dist(mlx, mly, tx, ty)) - get_bla() : -acos((tx-mlx)/dist(mlx, mly, tx, ty)) - get_bla();
}
double Dynein::get_ta() {
	return (mry > ty) ? acos((mrx-tx)/dist(tx, ty, mrx, mry)) - get_mla() - get_bla(): -acos((mrx-tx)/dist(tx, ty, mrx, mry)) - get_mla() - get_bla();
}
double Dynein::get_mra() {
	return (bry > mry) ? acos((brx-mrx)/dist(mrx, mry, brx, bry)) - get_ta() - get_mla() - get_bla() :
		-acos((brx-mrx)/dist(mrx, mry, brx, bry)) - get_ta() - get_mla() - get_bla();
}

void Dynein::log() {
	//Writing lFootX	lKneeX	hipX	rKneeX	rFootX	lFootY	lKneeY	hipY	rKneeY	rFootY	lHipAngle	rHipAngle	lKneeAngle	rKneeAngle
	
	FILE* file;
	file = fopen("data.txt", "w");
	
	fprintf(file, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t", blx, mlx, tx, mrx, brx, bly, mly, ty, mry, bry, get_ta(), get_mla(), get_mra());	
		
}

/*** Utility Functions ***/

double randAngle(double range) {
	srand(time(NULL));
	return ((((double)(rand() % 1000))/500) - 1)*range;
}

double dist(double x1, double y1, double x2, double y2) {
	return sqrt(pow(x2-x1,2)+pow(y2-y1,2));
}

Dynein* initProtein(Dynein* dyn) {
	dyn->set_blx(15);
	dyn->set_bly(15);
	
	dyn->set_bla(7.5*M_PI/18);
	dyn->set_mla(-4*M_PI/18);
	dyn->set_ta(-7*M_PI/18);
	dyn->set_mra(-4*M_PI/18);
}

int main() {
	Dynein* dyn = (Dynein*) malloc(sizeof(Dynein));
	initProtein(dyn);
	dyn->log();
	free(dyn);
	dyn = NULL;
	return 0;
}
