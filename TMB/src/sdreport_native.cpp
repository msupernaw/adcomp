#include <Eigen/Dense>
#include <R.h>
#include <Rinternals.h>

static SEXP get_list_element(SEXP list, const char *name)
{
    SEXP names = Rf_getAttrib(list, R_NamesSymbol);
    for (R_xlen_t index = 0; index < XLENGTH(list); ++index)
    {
        if (strcmp(CHAR(STRING_ELT(names, index)), name) == 0)
        {
            return VECTOR_ELT(list, index);
        }
    }
    return R_NilValue;
}

static Eigen::VectorXd as_vector(SEXP x)
{
    if (!Rf_isReal(x))
    {
        Rf_error("Expected numeric vector.");
    }
    const R_xlen_t size = XLENGTH(x);
    Eigen::VectorXd out(size);
    for (R_xlen_t index = 0; index < size; ++index)
    {
        out(index) = REAL(x)[index];
    }
    return out;
}

static Eigen::MatrixXd as_matrix(SEXP x)
{
    if (!Rf_isMatrix(x) || !Rf_isReal(x))
    {
        Rf_error("Expected numeric matrix.");
    }
    const R_xlen_t nrow = Rf_nrows(x);
    const R_xlen_t ncol = Rf_ncols(x);
    Eigen::MatrixXd out(nrow, ncol);
    for (R_xlen_t col = 0; col < ncol; ++col)
    {
        for (R_xlen_t row = 0; row < nrow; ++row)
        {
            out(row, col) = REAL(x)[row + nrow * col];
        }
    }
    return out;
}

static SEXP as_sexp(const Eigen::VectorXd &x)
{
    SEXP out = PROTECT(Rf_allocVector(REALSXP, x.size()));
    for (int index = 0; index < x.size(); ++index)
    {
        REAL(out)[index] = x(index);
    }
    UNPROTECT(1);
    return out;
}

static SEXP as_sexp(const Eigen::MatrixXd &x)
{
    SEXP out = PROTECT(Rf_allocMatrix(REALSXP, x.rows(), x.cols()));
    for (int col = 0; col < x.cols(); ++col)
    {
        for (int row = 0; row < x.rows(); ++row)
        {
            REAL(out)[row + x.rows() * col] = x(row, col);
        }
    }
    UNPROTECT(1);
    return out;
}

extern "C" SEXP tmb_sdreport_get_uncertainty(SEXP rep, SEXP name)
{
    if (!Rf_isNewList(rep))
    {
        Rf_error("'rep' must be an sdreport object.");
    }
    SEXP uncertainty = get_list_element(rep, "uncertainty");
    if (Rf_isNull(uncertainty))
    {
        Rf_error("No 'uncertainty' payload found. Run sdreport() on an object using SDREPORT_* registrations.");
    }
    if (Rf_isNull(name))
    {
        return uncertainty;
    }
    if (!Rf_isString(name) || Rf_length(name) != 1)
    {
        Rf_error("'name' must be NULL or a single string.");
    }
    const char *key = CHAR(STRING_ELT(name, 0));
    SEXP item = get_list_element(uncertainty, key);
    if (Rf_isNull(item))
    {
        Rf_error("Unknown uncertainty entry '%s'.", key);
    }
    return item;
}

extern "C" SEXP tmb_sdreport_get_scalar_estimate_sd(SEXP rep, SEXP name)
{
    if (!Rf_isString(name) || Rf_length(name) != 1)
    {
        Rf_error("'name' must be a single string.");
    }
    SEXP item = tmb_sdreport_get_uncertainty(rep, name);
    SEXP estimate = get_list_element(item, "estimate");
    SEXP sd = get_list_element(item, "sd");
    if (!Rf_isReal(estimate) || XLENGTH(estimate) != 1)
    {
        Rf_error("Entry '%s' estimate is not scalar.", CHAR(STRING_ELT(name, 0)));
    }
    if (!Rf_isReal(sd) || XLENGTH(sd) != 1)
    {
        Rf_error("Entry '%s' sd is not scalar.", CHAR(STRING_ELT(name, 0)));
    }
    SEXP out, out_names;
    PROTECT(out = Rf_allocVector(REALSXP, 2));
    PROTECT(out_names = Rf_allocVector(STRSXP, 2));
    REAL(out)[0] = REAL(estimate)[0];
    REAL(out)[1] = REAL(sd)[0];
    SET_STRING_ELT(out_names, 0, Rf_mkChar("estimate"));
    SET_STRING_ELT(out_names, 1, Rf_mkChar("sd"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    UNPROTECT(2);
    return out;
}

extern "C" SEXP tmb_sdreport_delta_fixed(SEXP phi_sexp, SEXP jac_sexp, SEXP Vtheta)
{
    if (!Rf_isReal(phi_sexp))
        Rf_error("'phi' must be numeric.");
    if (!Rf_isReal(jac_sexp) || !Rf_isMatrix(jac_sexp))
        Rf_error("'jac' must be a numeric matrix.");
    if (!Rf_isMatrix(Vtheta) || !Rf_isReal(Vtheta))
        Rf_error("'Vtheta' must be a numeric matrix.");

    Eigen::VectorXd phi = as_vector(phi_sexp);
    Eigen::MatrixXd jac = as_matrix(jac_sexp);
    Eigen::MatrixXd vtheta = as_matrix(Vtheta);

    if (jac.rows() != phi.size())
        Rf_error("Jacobian rows do not match phi length.");
    if (jac.cols() != vtheta.rows() || vtheta.rows() != vtheta.cols())
        Rf_error("Jacobian columns do not match Vtheta dimensions.");

    Eigen::MatrixXd cov = jac * vtheta * jac.transpose();
    Eigen::VectorXd sd(cov.rows());
    for (int i = 0; i < cov.rows(); ++i)
        sd(i) = sqrt(cov(i, i));

    SEXP ans, ans_names;
    PROTECT(ans = Rf_allocVector(VECSXP, 3));
    PROTECT(ans_names = Rf_allocVector(STRSXP, 3));
    SET_VECTOR_ELT(ans, 0, as_sexp(phi));
    SET_VECTOR_ELT(ans, 1, as_sexp(sd));
    SET_VECTOR_ELT(ans, 2, as_sexp(cov));
    SET_STRING_ELT(ans_names, 0, Rf_mkChar("value"));
    SET_STRING_ELT(ans_names, 1, Rf_mkChar("sd"));
    SET_STRING_ELT(ans_names, 2, Rf_mkChar("cov"));
    Rf_setAttrib(ans, R_NamesSymbol, ans_names);

    UNPROTECT(2);
    return ans;
}
