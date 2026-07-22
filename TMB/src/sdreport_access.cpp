#include "TMB.hpp"

extern "C" SEXP tmb_sdreport_get_uncertainty(SEXP rep, SEXP name)
{
    if (!Rf_isNewList(rep))
    {
        Rf_error("'rep' must be an sdreport object.");
    }
    SEXP uncertainty = getListElement(rep, "uncertainty");
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
    SEXP item = getListElement(uncertainty, key);
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
    SEXP estimate = getListElement(item, "estimate");
    SEXP sd = getListElement(item, "sd");
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
    REAL(out)
    [0] = REAL(estimate)[0];
    REAL(out)
    [1] = REAL(sd)[0];
    SET_STRING_ELT(out_names, 0, Rf_mkChar("estimate"));
    SET_STRING_ELT(out_names, 1, Rf_mkChar("sd"));
    Rf_setAttrib(out, R_NamesSymbol, out_names);
    UNPROTECT(2);
    return out;
}