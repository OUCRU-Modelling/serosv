// Generated by rstantools.  Do not edit by hand.

/*
    serosv is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    serosv is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with serosv.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_prevalence_correction_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 23> locations_array__ =
  {" (found before start of program)",
  " (in 'prevalence_correction', line 11, column 2 to column 31)",
  " (in 'prevalence_correction', line 12, column 2 to column 31)",
  " (in 'prevalence_correction', line 13, column 2 to column 36)",
  " (in 'prevalence_correction', line 16, column 9 to column 13)",
  " (in 'prevalence_correction', line 16, column 2 to column 30)",
  " (in 'prevalence_correction', line 18, column 2 to column 66)",
  " (in 'prevalence_correction', line 19, column 2 to column 66)",
  " (in 'prevalence_correction', line 22, column 6 to column 28)",
  " (in 'prevalence_correction', line 24, column 6 to column 68)",
  " (in 'prevalence_correction', line 26, column 6 to column 51)",
  " (in 'prevalence_correction', line 21, column 19 to line 28, column 3)",
  " (in 'prevalence_correction', line 21, column 2 to line 28, column 3)",
  " (in 'prevalence_correction', line 2, column 2 to column 20)",
  " (in 'prevalence_correction', line 3, column 2 to column 26)",
  " (in 'prevalence_correction', line 4, column 2 to column 26)",
  " (in 'prevalence_correction', line 5, column 2 to column 29)",
  " (in 'prevalence_correction', line 6, column 2 to column 29)",
  " (in 'prevalence_correction', line 7, column 20 to column 24)",
  " (in 'prevalence_correction', line 7, column 2 to column 26)",
  " (in 'prevalence_correction', line 8, column 18 to column 22)",
  " (in 'prevalence_correction', line 8, column 2 to column 24)",
  " (in 'prevalence_correction', line 13, column 30 to column 34)"};
#include <stan_meta_header.hpp>
class model_prevalence_correction final : public model_base_crtp<model_prevalence_correction> {
private:
  int Nage;
  double init_se;
  double init_sp;
  int study_size_se;
  int study_size_sp;
  std::vector<int> posi;
  std::vector<int> ni;
public:
  ~model_prevalence_correction() {}
  model_prevalence_correction(stan::io::var_context& context__, unsigned int
                              random_seed__ = 0, std::ostream*
                              pstream__ = nullptr) : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_prevalence_correction_namespace::model_prevalence_correction";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 13;
      context__.validate_dims("data initialization", "Nage", "int",
        std::vector<size_t>{});
      Nage = std::numeric_limits<int>::min();
      current_statement__ = 13;
      Nage = context__.vals_i("Nage")[(1 - 1)];
      current_statement__ = 13;
      stan::math::check_greater_or_equal(function__, "Nage", Nage, 1);
      current_statement__ = 14;
      context__.validate_dims("data initialization", "init_se", "double",
        std::vector<size_t>{});
      init_se = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 14;
      init_se = context__.vals_r("init_se")[(1 - 1)];
      current_statement__ = 14;
      stan::math::check_greater_or_equal(function__, "init_se", init_se, 0.5);
      current_statement__ = 15;
      context__.validate_dims("data initialization", "init_sp", "double",
        std::vector<size_t>{});
      init_sp = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 15;
      init_sp = context__.vals_r("init_sp")[(1 - 1)];
      current_statement__ = 15;
      stan::math::check_greater_or_equal(function__, "init_sp", init_sp, 0.5);
      current_statement__ = 16;
      context__.validate_dims("data initialization", "study_size_se", "int",
        std::vector<size_t>{});
      study_size_se = std::numeric_limits<int>::min();
      current_statement__ = 16;
      study_size_se = context__.vals_i("study_size_se")[(1 - 1)];
      current_statement__ = 16;
      stan::math::check_greater_or_equal(function__, "study_size_se",
        study_size_se, 0);
      current_statement__ = 17;
      context__.validate_dims("data initialization", "study_size_sp", "int",
        std::vector<size_t>{});
      study_size_sp = std::numeric_limits<int>::min();
      current_statement__ = 17;
      study_size_sp = context__.vals_i("study_size_sp")[(1 - 1)];
      current_statement__ = 17;
      stan::math::check_greater_or_equal(function__, "study_size_sp",
        study_size_sp, 0);
      current_statement__ = 18;
      stan::math::validate_non_negative_index("posi", "Nage", Nage);
      current_statement__ = 19;
      context__.validate_dims("data initialization", "posi", "int",
        std::vector<size_t>{static_cast<size_t>(Nage)});
      posi = std::vector<int>(Nage, std::numeric_limits<int>::min());
      current_statement__ = 19;
      posi = context__.vals_i("posi");
      current_statement__ = 19;
      stan::math::check_greater_or_equal(function__, "posi", posi, 0);
      current_statement__ = 20;
      stan::math::validate_non_negative_index("ni", "Nage", Nage);
      current_statement__ = 21;
      context__.validate_dims("data initialization", "ni", "int",
        std::vector<size_t>{static_cast<size_t>(Nage)});
      ni = std::vector<int>(Nage, std::numeric_limits<int>::min());
      current_statement__ = 21;
      ni = context__.vals_i("ni");
      current_statement__ = 21;
      stan::math::check_greater_or_equal(function__, "ni", ni, 0);
      current_statement__ = 22;
      stan::math::validate_non_negative_index("theta", "Nage", Nage);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = 1 + 1 + Nage;
  }
  inline std::string model_name() const final {
    return "model_prevalence_correction";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_prevalence_correction_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      local_scalar_t__ est_se = DUMMY_VAR__;
      current_statement__ = 1;
      est_se = in__.template read_constrain_lub<local_scalar_t__,
                 jacobian__>(0, 1, lp__);
      local_scalar_t__ est_sp = DUMMY_VAR__;
      current_statement__ = 2;
      est_sp = in__.template read_constrain_lub<local_scalar_t__,
                 jacobian__>(0, 1, lp__);
      std::vector<local_scalar_t__> theta =
        std::vector<local_scalar_t__>(Nage, DUMMY_VAR__);
      current_statement__ = 3;
      theta = in__.template read_constrain_lub<std::vector<local_scalar_t__>,
                jacobian__>(0, 1, lp__, Nage);
      {
        current_statement__ = 4;
        stan::math::validate_non_negative_index("apparent_theta", "Nage",
          Nage);
        Eigen::Matrix<local_scalar_t__,-1,1> apparent_theta =
          Eigen::Matrix<local_scalar_t__,-1,1>::Constant(Nage, DUMMY_VAR__);
        current_statement__ = 6;
        lp_accum__.add(stan::math::beta_lpdf<propto__>(est_se, (init_se *
                         study_size_se), ((1 - init_se) * study_size_se)));
        current_statement__ = 7;
        lp_accum__.add(stan::math::beta_lpdf<propto__>(est_sp, (init_sp *
                         study_size_sp), ((1 - init_sp) * study_size_sp)));
        current_statement__ = 12;
        for (int i = 1; i <= Nage; ++i) {
          current_statement__ = 8;
          lp_accum__.add(stan::math::beta_lpdf<propto__>(
                           stan::model::rvalue(theta, "theta",
                             stan::model::index_uni(i)), 1, 1));
          current_statement__ = 9;
          stan::model::assign(apparent_theta,
            ((stan::model::rvalue(theta, "theta", stan::model::index_uni(i))
            * est_se) + ((1 -
            stan::model::rvalue(theta, "theta", stan::model::index_uni(i))) *
            (1 - est_sp))), "assigning variable apparent_theta",
            stan::model::index_uni(i));
          current_statement__ = 10;
          lp_accum__.add(stan::math::binomial_lpmf<propto__>(
                           stan::model::rvalue(posi, "posi",
                             stan::model::index_uni(i)),
                           stan::model::rvalue(ni, "ni",
                             stan::model::index_uni(i)),
                           stan::model::rvalue(apparent_theta,
                             "apparent_theta", stan::model::index_uni(i))));
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_prevalence_correction_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      double est_se = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 1;
      est_se = in__.template read_constrain_lub<local_scalar_t__,
                 jacobian__>(0, 1, lp__);
      double est_sp = std::numeric_limits<double>::quiet_NaN();
      current_statement__ = 2;
      est_sp = in__.template read_constrain_lub<local_scalar_t__,
                 jacobian__>(0, 1, lp__);
      std::vector<double> theta =
        std::vector<double>(Nage, std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 3;
      theta = in__.template read_constrain_lub<std::vector<local_scalar_t__>,
                jacobian__>(0, 1, lp__, Nage);
      out__.write(est_se);
      out__.write(est_sp);
      out__.write(theta);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ est_se = DUMMY_VAR__;
      current_statement__ = 1;
      est_se = in__.read<local_scalar_t__>();
      out__.write_free_lub(0, 1, est_se);
      local_scalar_t__ est_sp = DUMMY_VAR__;
      current_statement__ = 2;
      est_sp = in__.read<local_scalar_t__>();
      out__.write_free_lub(0, 1, est_sp);
      std::vector<local_scalar_t__> theta =
        std::vector<local_scalar_t__>(Nage, DUMMY_VAR__);
      current_statement__ = 3;
      stan::model::assign(theta,
        in__.read<std::vector<local_scalar_t__>>(Nage),
        "assigning variable theta");
      out__.write_free_lub(0, 1, theta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "est_se", "double",
        std::vector<size_t>{});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "est_sp", "double",
        std::vector<size_t>{});
      current_statement__ = 3;
      context__.validate_dims("parameter initialization", "theta", "double",
        std::vector<size_t>{static_cast<size_t>(Nage)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      local_scalar_t__ est_se = DUMMY_VAR__;
      current_statement__ = 1;
      est_se = context__.vals_r("est_se")[(1 - 1)];
      out__.write_free_lub(0, 1, est_se);
      local_scalar_t__ est_sp = DUMMY_VAR__;
      current_statement__ = 2;
      est_sp = context__.vals_r("est_sp")[(1 - 1)];
      out__.write_free_lub(0, 1, est_sp);
      std::vector<local_scalar_t__> theta =
        std::vector<local_scalar_t__>(Nage, DUMMY_VAR__);
      current_statement__ = 3;
      theta = context__.vals_r("theta");
      out__.write_free_lub(0, 1, theta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"est_se", "est_sp", "theta"};
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {}
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{},
                std::vector<size_t>{},
                std::vector<size_t>{static_cast<size_t>(Nage)}};
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {}
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "est_se");
    param_names__.emplace_back(std::string() + "est_sp");
    for (int sym1__ = 1; sym1__ <= Nage; ++sym1__) {
      param_names__.emplace_back(std::string() + "theta" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {}
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    param_names__.emplace_back(std::string() + "est_se");
    param_names__.emplace_back(std::string() + "est_sp");
    for (int sym1__ = 1; sym1__ <= Nage; ++sym1__) {
      param_names__.emplace_back(std::string() + "theta" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {}
    if (emit_generated_quantities__) {}
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"est_se\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"est_sp\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nage) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"est_se\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"est_sp\",\"type\":{\"name\":\"real\"},\"block\":\"parameters\"},{\"name\":\"theta\",\"type\":{\"name\":\"array\",\"length\":" + std::to_string(Nage) + ",\"element_type\":{\"name\":\"real\"}},\"block\":\"parameters\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((1 + 1) + Nage);
    const size_t num_transformed = emit_transformed_parameters * (0);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((1 + 1) + Nage);
    const size_t num_transformed = emit_transformed_parameters * (0);
    const size_t num_gen_quantities = emit_generated_quantities * (0);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_prevalence_correction_namespace::model_prevalence_correction;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_prevalence_correction_namespace::profiles__;
}
#endif
#endif