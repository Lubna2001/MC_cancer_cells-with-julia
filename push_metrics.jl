# =========================================
# push_metrics.jl
# Counters to quantify pushing usage
# =========================================

module PushMetrics

export PushStats, reset!, record_div_success_space!, record_div_success_push!, record_push_attempt!, record_push_success!, metrics_summary

mutable struct PushStats
    # successful cancer divisions
    D_total::Int
    D_space::Int
    D_push::Int

    # pushing attempts and successes
    P_att::Int
    P_succ::Int

    function PushStats()
        new(0, 0, 0, 0, 0)
    end
end

function reset!(st::PushStats)
    st.D_total = 0
    st.D_space = 0
    st.D_push  = 0
    st.P_att   = 0
    st.P_succ  = 0
    return st
end

# Call when a division succeeded because an empty neighbor existed
function record_div_success_space!(st::PushStats)
    st.D_total += 1
    st.D_space += 1
    return nothing
end

# Call when a division succeeded after a push created space
function record_div_success_push!(st::PushStats)
    st.D_total += 1
    st.D_push  += 1
    return nothing
end

# Call every time you attempt a push
function record_push_attempt!(st::PushStats)
    st.P_att += 1
    return nothing
end

# Call when a push succeeds in creating space
function record_push_success!(st::PushStats)
    st.P_succ += 1
    return nothing
end

# Convenience summary numbers
function metrics_summary(st::PushStats)
    I_push = st.D_total == 0 ? 0.0 : st.D_push / st.D_total
    E_push = st.P_att   == 0 ? 0.0 : st.P_succ / st.P_att
    return (; D_total=st.D_total, D_space=st.D_space, D_push=st.D_push,
            P_att=st.P_att, P_succ=st.P_succ, I_push=I_push, E_push=E_push)
end

end # module
