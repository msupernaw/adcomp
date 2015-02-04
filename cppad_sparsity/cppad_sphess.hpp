namespace cppad_sparsity{
  using CppAD::vector;
  typedef vector< std::set<size_t> >  SetVector;
  typedef vector<bool>                BoolVector;
  /* Calculate pattern of sparse Jacobian */
  void calc_sparsity_jac(SetVector& sparsity_set, CppAD::ADFun<double>* pf)
  {	size_t n = pf->Domain();
    SetVector r_set(n);
    for(size_t j = 0; j < n; j++)
      r_set[j].insert(j);
    sparsity_set = pf->ForSparseJac(n, r_set);
  }
  /* Calculate pattern of sparse Hessian */
  void calc_sparsity_hes(SetVector& sparsity_set, CppAD::ADFun<double>* pf)
  {	size_t n = pf->Domain();
    size_t m = pf->Range();
    CPPAD_ASSERT_UNKNOWN( m == 1 );
    SetVector r_set(n);
    for(size_t i = 0; i < n; i++)
      r_set[i].insert(i);
    pf->ForSparseJac(n, r_set);
    //
    SetVector s_set(m);
    s_set[0].insert(0);
    //
    sparsity_set = pf->RevSparseHes(n, s_set);
  }
  /* Convert pattern to SEXP */
  SEXP asSEXP(std::set<size_t>& set){
    SEXP ans;
    PROTECT(ans = allocVector(INTSXP, set.size()));
    std::set<size_t>::iterator it;
    int k=0;
    for (it = set.begin(); it != set.end(); ++it){
      INTEGER(ans)[k++] = *it;
    }
    UNPROTECT(1);
    return ans;
  }
  SEXP asSEXP(SetVector& sparsity_set){
    SEXP ans;
    PROTECT(ans = allocVector(VECSXP, sparsity_set.size()));
    for(size_t j=0; j < sparsity_set.size(); j++){
      SET_VECTOR_ELT(ans, j, asSEXP(sparsity_set[j]));
    }
    UNPROTECT(1);
    return ans;
  }
}


extern "C"
{
  SEXP EvalADFunTest(SEXP f, SEXP theta, SEXP control)
  {
    if(!isNewList(control))error("'control' must be a list");
    typedef ADFun<double> ADFunType;  
    ADFunType* pf;
    pf=(ADFunType*)R_ExternalPtrAddr(f);
    PROTECT(theta=coerceVector(theta,REALSXP));
    int n=pf->Domain();
    if(LENGTH(theta)!=n)error("Wrong parameter length.");

    using namespace cppad_sparsity;
    SetVector s(n);
    calc_sparsity_hes(s, pf);

    UNPROTECT(1);
    return asSEXP(s);
  }

}
