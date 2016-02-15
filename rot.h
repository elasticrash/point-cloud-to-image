#include <math.h>

double r11(double phi, double kapa)
{
	return (cos(phi)*cos(kapa));
}
double r12(double omega, double phi, double kapa)
{
	return (cos(omega)*sin(kapa)+sin(omega)*sin(phi)*cos(kapa));
}
double r13(double omega, double phi, double kapa)
{
	return (sin(omega)*sin(kapa)-cos(omega)*sin(phi)*cos(kapa));
}
double r21(double phi, double kapa)
{
	return (-cos(phi)*sin(kapa));
}
double r22(double omega, double phi, double kapa)
{
	return (cos(omega)*cos(kapa)-sin(omega)*sin(phi)*sin(kapa));
}
double r23(double omega, double phi, double kapa)
{
	return (sin(omega)*cos(kapa)+cos(omega)*sin(phi)*sin(kapa));
}
double r31(double phi)
{
	return (sin(phi));
}
double r32(double omega, double phi)
{
	return (-sin(omega)*cos(phi));
}
double r33(double omega, double phi)
{
	return (cos(omega)*cos(phi));
}