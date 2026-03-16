function res = Proj(x,y,f)
    n = length(x); %#ok<NASGU>
    syms x_sym [n,1]
    grad = jacobian(f(x_sym),x_sym);
    grad_f = matlabFunction(grad, 'Vars', {x_sym});
    g = grad_f(x)';
    if y'*g > 0 && f(x) > 0
        res = y - g*g'/norm(g)^2*y*f(x);
    else
        res = y;
    end
end