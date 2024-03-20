#ifndef MATHFUNC_H
#define MATHFUNC_H

#define PI 3.14159265
double DeltaPhiSinFit(double DPhi, double a, double b, double c, double N) { return N * sin(a * DPhi - b) + c; }

double PowerLaw(double x, double a, double b, double c) { return a + b * pow(x, c); }

// the following two function ("Wild") are used to fit the parameters of the LevyTsallis. Dont ask about the name..
double WildFunction(double x, double a, double b, double c, double d, double e, double n, double m)
{
  return a * pow(x - m, 2) + b * x + c + d * pow(e * x, n);
}

double LevyTsallis(double x, double Norm, double n, double T, double m)
{
  return (x * Norm * (n - 1) * (n - 2)) / (n * T * (n * T + m * (n - 2))) *
         pow(1 + (sqrt(m * m + x * x) - m) / (n * T), -n);
}
double LevyTsallisTF1(double* x, double* params)
{
  double Norm = params[0];
  double n = params[1];
  double T = params[2];
  double m = params[3];
  return (x[0] * Norm * (n - 1) * (n - 2)) / (n * T * (n * T + m * (n - 2))) *
         pow(1 + (sqrt(m * m + x[0] * x[0]) - m) / (n * T), -n);
}

double GausSource(double r, double s)
{
  return 4 * PI * pow(r, 2) / pow(4 * PI * pow(s, 2), 3. / 2.) * exp(-pow(r, 2) / (4 * pow(s, 2)));
}
double NormalDist(double x, double mu, double sigma) { return exp(-0.5 * pow((x - mu) / (sigma), 2)); }

#endif