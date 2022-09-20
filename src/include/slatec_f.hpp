#ifndef SLATEC_F_H
#define SLATEC_F_H
//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

extern "C" {
// used by BspBasis class.
// creates the B-splines. deBoor routine
 void    dbspvd_(const double *t ,  const int &k ,
            const int &nderiv, const double &x ,
            const int &ileft,  const int &ldvnik,
            double *vnikx,     double *work);

double     bvalue_(double * t, double * c, const int & n, 
            const int & k, const double & x, const int & d) ;

double     bvalu_(double * t, double * c, const int & n, 
            const int & k, const int & p, const double & x, const int & d, double * work) ;
}

#endif // SLATEC_F_H