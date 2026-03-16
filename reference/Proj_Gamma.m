function res = Proj_Gamma(Gamma,x,y,f)
    n = length(x); %#ok<NASGU>
    syms x_sym [n,1]
    grad = jacobian(f(x_sym),x_sym);
    grad_f = matlabFunction(grad, 'Vars', {x_sym});
    g = grad_f(x)';
    if y'*g > 0 && f(x) > 0
        res = Gamma*y - Gamma*g*g'/(g'*Gamma*g)*Gamma*y*f(x);
    else
        res = Gamma*y;
    end
end