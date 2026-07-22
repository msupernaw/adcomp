#ifndef TMB_SDREPORT_ACCESS_HPP
#define TMB_SDREPORT_ACCESS_HPP

#include <R_ext/Rdynload.h>
#include <map>
#include <string>
#include <utility>
#include <vector>

struct sdreport_array
{
    std::vector<int> dim;
    std::vector<double> value;

    bool is_scalar() const
    {
        return value.size() == 1;
    }

    double scalar() const
    {
        if (!is_scalar())
        {
            Rf_error("Requested scalar from non-scalar sdreport value.");
        }
        return value[0];
    }

    vector<double> as_vector() const
    {
        vector<double> out(value.size());
        for (size_t i = 0; i < value.size(); ++i)
            out[i] = value[i];
        return out;
    }

    matrix<double> as_matrix() const
    {
        if (dim.size() != 2)
        {
            Rf_error("Requested matrix from sdreport value with rank %d.", (int)dim.size());
        }
        vector<double> out = as_vector();
        return asMatrix(out, dim[0], dim[1]);
    }
};

inline std::vector<int> sdreport_dim(SEXP x)
{
    SEXP dim = Rf_getAttrib(x, R_DimSymbol);
    std::vector<int> out;
    if (Rf_isNull(dim))
    {
        out.push_back(Rf_length(x));
        return out;
    }
    out.resize(Rf_length(dim));
    for (int i = 0; i < Rf_length(dim); ++i)
        out[i] = INTEGER(dim)[i];
    return out;
}

inline std::vector<double> sdreport_values(SEXP x)
{
    if (!Rf_isReal(x))
    {
        Rf_error("Expected numeric sdreport value.");
    }
    R_xlen_t n = XLENGTH(x);
    std::vector<double> out(n);
    double *ptr = REAL(x);
    for (R_xlen_t i = 0; i < n; ++i)
        out[i] = ptr[i];
    return out;
}

struct sdreport_entry
{
    std::string name;
    std::string type;
    sdreport_array estimate;
    sdreport_array std_error;
};

struct sdreport_access
{
    std::vector<sdreport_entry> entries;
    std::map<std::string, size_t> index;

    sdreport_access() {}

    explicit sdreport_access(SEXP x)
    {
        assign(x);
    }

    void clear()
    {
        entries.clear();
        index.clear();
    }

    void assign(SEXP x)
    {
        clear();
        SEXP uncertainty = getListElement(x, "uncertainty");
        if (Rf_isNull(uncertainty))
            return;
        SEXP names = Rf_getAttrib(uncertainty, R_NamesSymbol);
        for (int i = 0; i < Rf_length(uncertainty); ++i)
        {
            SEXP item = VECTOR_ELT(uncertainty, i);
            sdreport_entry entry;
            entry.name = CHAR(STRING_ELT(getListElement(item, "name"), 0));
            entry.type = CHAR(STRING_ELT(getListElement(item, "type"), 0));
            SEXP estimate = getListElement(item, "estimate");
            SEXP std_error = getListElement(item, "sd");
            entry.estimate.dim = sdreport_dim(estimate);
            entry.estimate.value = sdreport_values(estimate);
            entry.std_error.dim = sdreport_dim(std_error);
            entry.std_error.value = sdreport_values(std_error);
            entries.push_back(entry);
            std::string key = CHAR(STRING_ELT(names, i));
            index[key] = entries.size() - 1;
        }
    }

    bool has(const std::string &name) const
    {
        return index.find(name) != index.end();
    }

    const sdreport_entry &get(const std::string &name) const
    {
        std::map<std::string, size_t>::const_iterator it = index.find(name);
        if (it == index.end())
        {
            Rf_error("Unknown sdreport uncertainty entry '%s'.", name.c_str());
        }
        return entries[it->second];
    }
};

typedef SEXP (*tmb_sdreport_get_uncertainty_t)(SEXP, SEXP);
typedef SEXP (*tmb_sdreport_get_scalar_estimate_sd_t)(SEXP, SEXP);

inline tmb_sdreport_get_uncertainty_t tmb_sdreport_getter()
{
    static tmb_sdreport_get_uncertainty_t fn = NULL;
    if (fn == NULL)
    {
        fn = (tmb_sdreport_get_uncertainty_t)
            R_GetCCallable("TMB", "tmb_sdreport_get_uncertainty");
    }
    return fn;
}

inline tmb_sdreport_get_scalar_estimate_sd_t tmb_sdreport_scalar_getter()
{
    static tmb_sdreport_get_scalar_estimate_sd_t fn = NULL;
    if (fn == NULL)
    {
        fn = (tmb_sdreport_get_scalar_estimate_sd_t)
            R_GetCCallable("TMB", "tmb_sdreport_get_scalar_estimate_sd");
    }
    return fn;
}

inline SEXP tmb_get_sdreport_uncertainty(SEXP rep, const char *name = NULL)
{
    if (name == NULL)
    {
        return tmb_sdreport_getter()(rep, R_NilValue);
    }
    SEXP key;
    PROTECT(key = Rf_mkString(name));
    SEXP ans = tmb_sdreport_getter()(rep, key);
    UNPROTECT(1);
    return ans;
}

inline std::pair<double, double> tmb_get_sdreport_scalar(SEXP rep, const char *name)
{
    if (name == NULL)
    {
        Rf_error("'name' must not be NULL for scalar sdreport lookup.");
    }
    SEXP key;
    PROTECT(key = Rf_mkString(name));
    SEXP ans = tmb_sdreport_scalar_getter()(rep, key);
    if (!Rf_isReal(ans) || XLENGTH(ans) != 2)
    {
        UNPROTECT(1);
        Rf_error("Unexpected return from scalar sdreport accessor.");
    }
    std::pair<double, double> out(REAL(ans)[0], REAL(ans)[1]);
    UNPROTECT(1);
    return out;
}

#endif