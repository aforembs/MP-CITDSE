#ifndef MODEL_V_H
#define MODEL_V_H

#include <cassert>
#include <vector>
#include <au.h>
#include <string>

//xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//Model Potentials and parameters taken from Notre Dame Group:
//W.R.Johnson,D.S.Guo,M.Idrees and J.Sapisrtein
// Phys.Rev A(34),1043,1986   and Phys.Rev A(32),2093,1985
// University of Notre Dame ,Notre Dame ,Indiana
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


inline  double V_Rel_1_r( double x ){
  //
  assert(x) ;
  double aNucl ;
	 //	 double Vnucl ;
  aNucl = 2.2677e-5 ;
  
  //	 Vnucl = 0.5 *( (x*x) / (aNucl*aNucl)-3) /aNucl;

  if( x<aNucl )        return ( 1./x ) ;         //  (-Vnucl);
  else 		       return  ( 1./x ) ;
} 

// Abstract Base Class ( ABC ) for model potentials

class ModelV {
	 
public:
  virtual double V(const double & x) = 0 ;
  virtual void   plotV() = 0 ;
  virtual double V_InsideNucleus(const double& x, const double& aNucl) {
	  return  0.5 * ( (x*x) / (aNucl*aNucl)-3)/aNucl ;
	} 
	virtual ~ModelV() {} 
} ;

class V_c: public ModelV {
	private:
		const double c_ ;

	public :
		V_c(const double& c) : c_(c) {} 

		virtual double V([[maybe_unused]] const double& x)  { return  c_ ; } 
		virtual void   plotV() {}
};

class V_1_r: public ModelV {
	private:
		const double  z_ ;              // atomic number  
	 
	public :
		V_1_r(const double& z) : z_(z) {} 

		virtual inline double V(const double& x)  {       	 
		  return  ( z_) / x  ;
		} 
	 	virtual void   plotV() {}
} ;

class V_c_r2: public ModelV {
	private:
		const double c_ ;

	public :
		V_c_r2(const double& c) : c_(c) {} 

		virtual inline double V(const double& x)	{      	 
		  return   c_ / (x*x)  ;
	 	} 
		virtual void   plotV(){}
};

class V_r: public ModelV {
	public :
		V_r() {} 

		virtual inline double V(const double& x)  {  return  x ;}
		virtual void   plotV(){}
};

class V_r2: public ModelV {
	public :
		V_r2() {} 

		virtual inline double V(const double& x)  {  return  x*x ;}
		virtual void   plotV(){}
};

#endif


