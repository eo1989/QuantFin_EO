## Finding the implied volatility w/ Newton's Method

using ForwardDiff, Distributions

## Setting up the BSM formulae

d = Normal()
d1(S, K, T, r, v) = (log(S / K) + (r + v * v / 2) * T) / (v * √(T))
d2(S, K, T, r, v) = (log(S / K) + (r + v * v / 2) * T) / (v * √(T)) - v * √(T)

call_price(S, K, T, r, v, q) = S * exp(-q * T) * cdf(d, d1(S, K, T, r, v)) - K * exp(-r * T) * cdf(d, d2(S, K, T, r, v))

put_price(S, K, T, r, v, q) = K * exp(-r * T) * cdf(d, -d2(S, K, T, r, v)) - S * exp(-q * T) * cdf(d, -d1(S, K, T, r, v))

## Newton-Raphson method for finding roots
# could just use the Roots library..

function iter_newton(𝞼₀, tolerance, max_iter, mkt_price, S, K, T, r, q=0)
  # 𝞼₀: initial guess for σ, 0.5 recommended
  # tolerance: precision of the approximation
  # max_iter: set max iteration of the algo
  # mkt_price: market price of the option
  # S: spot price
  # K: strike price
  # T: time to maturity
  # r: interest rate
  # q: dividend paid, set to 0

  # Intrinsic value of dITM options may prevent the algorithm from converging.
  # May need to use put-call parity & find the volatility of the (opposite) OTM option instead.

    if K >= S
        price_v = v -> call_price(S, K, T, r, v, q)
    else
        r = -r
        price_v = v -> put_price(S, K, T, r, v, q)
    end

    Vega(v) = ForwardDiff.derivative(price_v, v)
    iter = 0
    v = 𝞼₀
    𝜟 = price_v(v) - mkt_price
    while abs(𝜟) > tolerance && iter < max_iter
        v = v - 0.1(price_v(v) - mkt_price / Vega(v))
        𝜟 = price_v(v) - mkt_price
        iter = iter + 1
        print("volatility was $v") 
    end
    return v, iter
end

#  Ex: =#
#  S = 200 =#
#  K = 250 =#
#  T = 97 =#
#  r = 0.0134 =#
#  q = 0 =#
#  @time iter_newton(0.5, 0.000001, 1000, 0.01, 183, 1.4, 0.0379, 0)
