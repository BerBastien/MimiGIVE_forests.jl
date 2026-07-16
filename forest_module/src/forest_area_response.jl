using Mimi

# =============================================================================
# Component B:  ForestAreaResponse
# -----------------------------------------------------------------------------
# Converts GIVE global-mean surface temperature (GMST) into an impact-region
# forest-area projection, relative to a per-run 2022 baseline.
#
#   delta_T_t = GMST_t - GMST_2022        (re-centred on THIS run's own 2022 value)
#   change_{r,t} = a_r + b1_r*delta_T_t + b2_r*delta_T_t^2   (units set by config)
#
#   :percent      A_{r,t} = A_{r,2022} * (1 + change/100)
#   :proportion   A_{r,t} = A_{r,2022} * (1 + change)
#   :log_change   A_{r,t} = A_{r,2022} * exp(change)
#
# Physical constraints:  A >= 0  (lower clip), and optionally A <= total region
# area (upper clip).  The raw, unconstrained projection is retained for diagnostics.
#
# IMPORTANT (marginal model): because this component runs *inside* each GIVE model
# instance, the base model and the pulse model each re-derive their OWN GMST_2022
# from their OWN FaIR trajectory. That is exactly what the spec requires --- the
# 2022 reference is never taken from an external observational number.
#
# R analogy: `mutate(dT = gmst - gmst[year==2022])` computed separately within each
# scenario's data frame, then a vectorised polynomial applied per region.
# =============================================================================

@defcomp ForestAreaResponse begin

    impact_region = Index()

    # ---- Inputs -------------------------------------------------------------
    gmst = Parameter(index=[time], unit="degC")   # connect to GIVE :temperature => :T

    beta1     = Parameter(index=[impact_region])            # beta_delta_gmst
    beta2     = Parameter(index=[impact_region])            # beta_delta_gmst_squared
    intercept = Parameter(index=[impact_region])            # alpha_r (0 if none supplied)

    baseline_forest_area     = Parameter(index=[impact_region])   # A_{r,2022}, in FOREST_AREA_UNITS
    total_impact_region_area = Parameter(index=[impact_region])   # used only if apply_upper_clip

    # Optional fitted delta-T range per region for extrapolation flagging.
    # Defaults (set at build) of -Inf / +Inf disable flagging.
    delta_gmst_fit_min = Parameter(index=[impact_region])
    delta_gmst_fit_max = Parameter(index=[impact_region])

    # ---- Configuration (scalars) -------------------------------------------
    baseline_year     = Parameter{Int}(default=2022)
    change_units_code = Parameter{Int}(default=1)   # 1=percent, 2=proportion, 3=log_change
    apply_lower_clip  = Parameter{Bool}(default=true)
    apply_upper_clip  = Parameter{Bool}(default=false)

    # ---- Outputs ------------------------------------------------------------
    gmst_baseline           = Variable()                              # GMST in the baseline year (scalar, this run)
    delta_gmst              = Variable(index=[time], unit="degC")
    delta_gmst_squared      = Variable(index=[time], unit="degC^2")
    forest_change           = Variable(index=[time, impact_region])   # raw fitted change, in coefficient units
    projected_forest_area_raw = Variable(index=[time, impact_region]) # unconstrained
    projected_forest_area   = Variable(index=[time, impact_region])   # constrained (>=0, optional cap)
    forest_area_clipped_lower = Variable(index=[time, impact_region]) # 1 if clipped at 0
    forest_area_clipped_upper = Variable(index=[time, impact_region]) # 1 if clipped at total area
    extrapolation_flag        = Variable(index=[time, impact_region]) # 1 if delta_T outside fitted range

    function init(p, v, d)
        # NaN sentinel: if the baseline year is never reached the deltas become
        # NaN (a loud failure) rather than silently wrong. Build-time validation
        # guarantees the baseline year is within the component's run range.
        v.gmst_baseline = NaN
    end

    function run_timestep(p, v, d, t)
        year = gettime(t)

        # Capture this run's own GMST in the baseline year.
        if year == p.baseline_year
            v.gmst_baseline = p.gmst[t]
        end

        # delta_T is defined from the baseline year onward; 0 before it.
        dT = year >= p.baseline_year ? (p.gmst[t] - v.gmst_baseline) : 0.0
        v.delta_gmst[t]         = dT
        v.delta_gmst_squared[t] = dT^2

        for ir in d.impact_region
            change = p.intercept[ir] + p.beta1[ir] * dT + p.beta2[ir] * dT^2
            v.forest_change[t, ir] = change

            A0 = p.baseline_forest_area[ir]
            raw = p.change_units_code == 1 ? A0 * (1 + change / 100) :
                  p.change_units_code == 2 ? A0 * (1 + change)       :
                  p.change_units_code == 3 ? A0 * exp(change)        :
                  error("Unknown change_units_code = $(p.change_units_code)")
            v.projected_forest_area_raw[t, ir] = raw

            A = raw
            lower = 0
            upper = 0
            if p.apply_lower_clip && A < 0
                A = 0.0
                lower = 1
            end
            if p.apply_upper_clip
                cap = p.total_impact_region_area[ir]
                if cap > 0 && A > cap
                    A = cap
                    upper = 1
                end
            end
            v.projected_forest_area[t, ir]     = A
            v.forest_area_clipped_lower[t, ir] = lower
            v.forest_area_clipped_upper[t, ir] = upper

            v.extrapolation_flag[t, ir] =
                (dT < p.delta_gmst_fit_min[ir] || dT > p.delta_gmst_fit_max[ir]) ? 1 : 0
        end
    end
end
