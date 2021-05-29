#include "CG/CGSolver.h"

CGSolver::CGSolver():
  is_initialized(false),
  toler(1.e-12)
{
  dealii::SparseMatrix<double> Temp;
  A=&Temp;
  //res.reinit(0);
  //Ap.reinit(0);
};

CGSolver::CGSolver(const Matrix& M, double tol)
{
  toler=tol;
  A=&M;
};

void CGSolver::reinit(const Matrix& A)
{
  is_initialized=true;
};

CGSolver::~CGSolver()
{
  clear();
}

void CGSolver::clear()
{
  
}

void CGSolver::get_Ap(std::vector<double>x)
{
  // if Ap is std::vector, then this following two commands are useful;
  //Ap.clear();
  //Ap.resize(x.size());
  Ap.reinit(x.size());
  if (A->n()!=Ap.size())
  {
    std::cout<<"Function multiply() ERROR: The sizes of matrix and vector are not same! Please check it!\n";
  }
  else
  {
    //std::cout<<"This is CGSolver::get_Ap!!! \n";
    for(int k=0;k<A->m();k++)
      {
	//	std::cout<<"this is the "<<k<<"th iterations \n";
	dealii::SparseMatrix<double>::const_iterator i=A->begin(k);
	Ap[k]=0;
        while(i!=A->end(k))
	  {
	    Ap[k]+=i->value()*x[i->column()];
	    ++i;
	  }
      }
  }
}

// This is the function for Ap function for dealii::Vector
void CGSolver::get_Ap(dealii::Vector<double>x)
{
  Ap.reinit(x.size());
  if (A->n()!=Ap.size())
  {
    std::cout<<"Function multiply() ERROR: The sizes of matrix and vector are not same! Please check it!\n";
  }
  else
  {
    for(int k=0;k<A->m();k++)
      {
	dealii::SparseMatrix<double>::const_iterator i=A->begin(k);
	Ap[k]=0;
        while(i!=A->end(k))
	  {
	    Ap[k]+=i->value()*x[i->column()];
	    ++i;
	  }
      }
  }
}

void CGSolver::get_res(std::vector<double> x, std::vector<double> r)
{
  res.reinit(x.size());
  for (int k=0;k<res.size();k++)
    {
      res[k]=0;
      dealii::SparseMatrix<double>::const_iterator i=A->begin(k);
      while (i!=A->end(k))
	{
	  res[k]+=i->value()*x[i->column()];
	  ++i;
	}
      res[k]=r[k]-res[k];
    }
}

// This is get_res function for dealii::Vector
void CGSolver::get_res(const dealii::Vector<double> x, const dealii::Vector<double> r)
{
  res.reinit(x.size());
  
  for (int k=0;k<res.size();k++)
    {
      res[k]=0;
      dealii::SparseMatrix<double>::const_iterator i=A->begin(k);
      while (i!=A->end(k))
	{
	  res[k]+=i->value()*x[i->column()];
	  ++i;
	}
      res[k]=r[k]-res[k];
    }
}


// frorbenius_norm function for 1-dimension vector:
// i.e. || a || _ F;
double frobenius_norm(std::vector<double> a)
{
  double norm=0;
  for(int k=0;k<a.size();k++)
    {
      norm+=a[k]*a[k];
    }
  return sqrt(norm);
}

// inner_product of two 1-dimension vector:
// i.e. <a, b>;
double inner_product(std::vector<double> a, std::vector<double> b)
{
  double sum=0;
  assert(a.size()==b.size());
  for(int k=0; k<a.size(); k++)
    {
      sum+=a[k]*b[k];
    }
  return sum;
}
// inner_product of two 1-dimension vector:
// i.e. <a, b>;
double inner_product(dealii::Vector<double> a, dealii::Vector<double> b)
{
  double sum=0;
  assert(a.size()==b.size());
  for(int k=0; k<a.size(); k++)
    {
      sum+=a[k]*b[k];
    }
  return sum;
}
// inner_product of two 1-dimension vector:
// i.e. <a, b>;
double inner_product(std::vector<double> a, dealii::Vector<double> b)
{
  double sum=0;
  assert(a.size()==b.size());
  for(int k=0; k<a.size(); k++)
    {
      sum+=a[k]*b[k];
    }
  return sum;
}

// max_norm:
double max_norm(std::vector<double> x)
{
  double max=0;
  for(int i=0;i<x.size();i++)
    {
      if(max<x[i])
	{
	  max=x[i];
	}
    }
  return max;
}
// max_norm:
double max_norm(dealii::Vector<double> x)
{
  double max=0;
  for(int i=0;i<x.size();i++)
    {
      if(max<fabs(x[i]))
	{
	  max=fabs(x[i]);
	}
    }
  return max;
}


void CGSolver::solve(std::vector<double>& x, const std::vector<double> r, double tol, int max_iter)
{
  Assert(is_initialized == true, ExcNotInitialized());
  if(tol==0.0)tol=toler;

  get_res(x,r);

  std::vector<double> temp_res(res.size());
  for(int i=0;i<res.size();i++)
    {
      temp_res[i]=res[i];
    }
  get_Ap(temp_res);
  
  double beta=0;
  double delta=0;

  std::vector<double> p(x.size(),0);
  for(int k=0;k<p.size();k++)
    {
      p[k]=-res[k];
    }

  get_Ap(p);
  
  delta=inner_product(p,res)/inner_product(p,Ap);
  // update vector x for the first iteration;
  for(int k=0;k<x.size();k++)
    {
      x[k]=x[k]+delta*p[k];
    }
  
  int iter=1;
  get_res(x,r);

  // in this iteration the extinction condition;
  // the max_norm can be replaced by frobenius_norm;
  while(iter<=max_iter&& max_norm(res)>tol)
    {
      //get_res(x,r);
      beta=inner_product(res,Ap)/inner_product(p,Ap);
      //beta=-beta;
      for(int k=0;k<p.size();k++)
	{
	  p[k]=-res[k]+beta*p[k];
	}

      //update Ap;
      get_Ap(p);

      delta=inner_product(p,res)/inner_product(p,Ap);

      // update vector x for the first iteration;
      for(int k=0;k<x.size();k++)
	{
	  x[k]=x[k]+delta*p[k];
	}
      get_res(x,r);
      iter++;
    }
 
}

void CGSolver::solve(dealii::Vector<double>& x, const dealii::Vector<double>& r, double tol, int max_iter)
{
  std::cout<<" This is CGSolver for dealii !!!!\n";
  
  Assert(is_initialized == true, ExcNotInitialized());
  if(tol==0.0)tol=toler;
  get_res(x,r);
  // p0=-g;
  get_Ap(res);
  double beta=0;
  double delta=0;
  //beta=inner_product(res,Ap)/inner_product(res,Ap);
  //beta=-beta;
  //std::vector<double> p(x.size(),0);
  dealii::Vector<double> p(x.size());
  for(int k=0;k<p.size();k++)
    {
      p[k]=-res[k];
    }
  get_Ap(p);
  delta=inner_product(p,res)/inner_product(p,Ap);

  // update vector x for the first iteration;
  for(int k=0;k<x.size();k++)
    {
      x[k]=x[k]+delta*p[k];
    }
  get_res(x,r);
  int iter=1;

  // in this iteration the extinction condition;
  // the max_norm can be replaced by frobenius_norm;
  while(iter<=max_iter&& max_norm(res)>tol)
    {

      beta=inner_product(res,Ap)/inner_product(p,Ap);
      //beta=-beta;
      for(int k=0;k<p.size();k++)
	{
	  p[k]=-res[k]+beta*p[k];
	}

      //update Ap;
      get_Ap(p);

      delta=inner_product(p,res)/inner_product(p,Ap);

      // update vector x for the first iteration;
      for(int k=0;k<x.size();k++)
	{
	  x[k]=x[k]+delta*p[k];
	}
      get_res(x,r);

      iter++;
    }

  std::cerr << "CGSolver converge at "<<iter-1<<" step"<<std::endl;
 
}
