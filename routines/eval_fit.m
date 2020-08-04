function [tret, vmax, erosion, totvol] = eval_fit(fitres, v_ret, t_ret, v_min, v_max)

a = fitres.a;
b = fitres.b;

tret = 1/(a) * v_ret^(-b);
vmax = (t_ret * a)^(-1.0/b);

if (v_max < 0) 
    v_max = vmax;
end


totvol = -a*b/(1+b) * (v_max^(1+b) - v_min^(1+b));
erosion = totvol/5.2e6;
