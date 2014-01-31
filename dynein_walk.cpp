#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <math.h>

const double thighLen = 10;
const double calfLen = 10;

class Dynein {
	public:
		void set_hipX(double d), set_hipY(double d);
		void set_lKneeX(double d), set_lKneeY(double d);
		void set_rKneeX(double d), set_rKneeY(double d);
		void set_lFootX(double d), set_lFootY(double d);
		void set_rFootX(double d), set_rFootY(double d);
		
		void set_lHipAngle(double d), set_rHipAngle(double d);
		void set_lKneeAngle(double d), set_rKneeAngle(double d);
		
		double get_hipX(), get_hipY();
		double get_lKneeX(), get_lKneeY();
		double get_rKneeX(), get_rKneeY();
		double get_lFootX(), get_lFootY();
		double get_rFootX(), get_rFootY();
		
		double get_lHipAngle(), get_rHipAngle();
		double get_lKneeAngle(), get_rKneeAngle();
		
	private:
		double hipX, hipY;
		double lKneeX, lKneeY;
		double rKneeX, rKneeY;
		double lFootX, lFootY;
		double rFootX, rFootY;
};

/*** Accessor Methods ***/

void Dynein::set_hipX(double d) {
	hipX = d;
}

void Dynein::set_hipY(double d) {
	hipY = d;
}

void Dynein::set_lKneeX(double d) {
	lKneeX = d;
}

void Dynein::set_rKneeX(double d) {
	rKneeX = d;
}

void Dynein::set_lKneeY(double d) {
	lKneeY = d;
}

void Dynein::set_rKneeY(double d) {
	rKneeY = d;
}

void Dynein::set_lFootX(double d) {
	lFootX = d;
}

void Dynein::set_rFootX(double d) {
	rFootX = d;
}

void Dynein::set_lFootY(double d) {
	lFootY = d;
}

void Dynein::set_rFootY(double d) {
	rFootY = d;
}

void Dynein::set_lHipAngle(double d) {
	double dX = hipX + thighLen*cos(d) - lKneeX;
	double dY = hipY + thighLen*sin(d) - lKneeY;
	lKneeX = hipX + thighLen*cos(d);
	lKneeY = hipY + thighLen*sin(d);
	lFootX += dX;
	lFootY += dY;
}

void Dynein::set_rHipAngle(double d) {
	double dX = hipX + thighLen*cos(d) - rKneeX;
	double dY = hipY + thighLen*sin(d) - rKneeY;
	rKneeX = hipX + thighLen*cos(d);
	rKneeY = hipY + thighLen*sin(d);
	rFootX += dX;
	rFootY += dY;
}

void Dynein::set_lKneeAngle(double d) {
	lFootX = lKneeX + calfLen*cos(d);
	lFootY = lKneeY + calfLen*sin(d);
}

void Dynein::set_rKneeAngle(double d) {
	rFootX = rKneeX + calfLen*cos(d);
	rFootY = rKneeY + calfLen*sin(d);
}

double Dynein::get_hipX() {
	return hipX;
}

double Dynein::get_hipY() {
	return hipY;
}

double Dynein::get_lKneeY() {
	return lKneeY;
}

double Dynein::get_rKneeY() {
	return rKneeY;
}

double Dynein::get_lKneeX() {
	return lKneeX;
}

double Dynein::get_rKneeX() {
	return rKneeX;
}

double Dynein::get_lFootY() {
	return lFootY;
}

double Dynein::get_rFootY() {
	return rFootY;
}

double Dynein::get_lFootX() {
	return lFootX;
}

double Dynein::get_rFootX() {
	return rFootX;
}

double Dynein::get_lHipAngle() {
	double hyp = sqrt(pow((hipX - lKneeX),2) + pow((hipY - lKneeY),2));
	double adj = lKneeX - hipX;
	if (hipY <= lKneeY) return acos(adj/hyp);
	else return -acos(adj/hyp);
}

double Dynein::get_rHipAngle() {
	double hyp = sqrt(pow((hipX - rKneeX),2) + pow((hipY - rKneeY),2));
	double adj = rKneeX - hipX;
	if (hipY <= rKneeY) return acos(adj/hyp);
	else return -acos(adj/hyp);
}

double Dynein::get_lKneeAngle() {
	double hyp = sqrt(pow((lKneeX - lFootX),2) + pow((lKneeY - lFootY),2));
	double adj = lFootX - lKneeX;
	if (lKneeY <= lFootY) return acos(adj/hyp);
	else return -acos(adj/hyp);
}

double Dynein::get_rKneeAngle() {
	double hyp = sqrt(pow((rKneeX - rFootX),2) + pow((rKneeY - rFootY),2));
	double adj = rFootX - rKneeX;
	if (rKneeY <= rFootY) return acos(adj/hyp);
	else return -acos(adj/hyp);
}

/*** Utility Functions ***/

double randAngle(double range) {
	srand(time(NULL));
	return ((((double)(rand() % 1000))/500) - 1)*range;
}

Dynein* initProtein(Dynein* dyn) {
	
	srand(time(NULL));
	
	dyn->set_hipX(0);
	dyn->set_hipY(20);
	
	dyn->set_lHipAngle(-5*M_PI/6 + randAngle(1));
	dyn->set_rHipAngle(-M_PI/6 + randAngle(1));
	
	dyn->set_lKneeAngle(-M_PI/3 + randAngle(1));
	dyn->set_rKneeAngle(-2*M_PI/3 + randAngle(1));
}

void printProtein(Dynein* dyn) {
	printf("hipX \t %f\n", dyn->get_hipX());
	printf("hipY \t %f\n", dyn->get_hipY());
	
	printf("lKneeX \t %f\n", dyn->get_lKneeX());
	printf("lKneeY \t %f\n", dyn->get_lKneeY());
	
	printf("rKneeX \t %f\n", dyn->get_rKneeX());
	printf("rKneeY \t %f\n", dyn->get_rKneeY());
	
	printf("lFootX \t %f\n", dyn->get_lFootX());
	printf("lFootY \t %f\n", dyn->get_lFootY());
	
	printf("rFootX \t %f\n", dyn->get_rFootX());
	printf("rFootY \t %f\n", dyn->get_rFootY());
	
	printf("lHipAngle: \t %f\n", dyn->get_lHipAngle());
	printf("rHipAngle: \t %f\n", dyn->get_rHipAngle());
	
	printf("lKneeAngle: \t %f\n", dyn->get_lKneeAngle());
	printf("rKneeAngle: \t %f\n", dyn->get_rKneeAngle());
}

void logProtein(Dynein* dyn) {
	//Writing lFootX	lKneeX	hipX	rKneeX	rFootX	lFootY	lKneeY	hipY	rKneeY	rFootY	lHipAngle	rHipAngle	lKneeAngle	rKneeAngle
	
	FILE* file;
	file = fopen("data.txt", "w");
	
	fprintf(file, "%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f", dyn->get_lFootX(), dyn->get_lKneeX(), dyn->get_hipX(),
		dyn->get_rKneeX(), dyn->get_rFootX(), dyn->get_lFootY(), dyn->get_lKneeY(), dyn->get_hipY(), dyn->get_rKneeY(),
		dyn->get_rFootY(), dyn->get_lHipAngle(), dyn->get_rHipAngle(), dyn->get_lKneeAngle(), dyn->get_rKneeAngle());	
		
}

int main() {
	Dynein* dyn = (Dynein*) malloc(sizeof(Dynein));
	initProtein(dyn);
	logProtein(dyn);
	free(dyn);
	dyn = NULL;
	return 0;
}
