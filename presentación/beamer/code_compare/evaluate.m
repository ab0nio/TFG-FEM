function total = evaluate(n)
% Evalúa la función de base xi*eta^2*zeta^2 sobre una nube n x n x n
% de puntos del hexaedro de referencia [0,1]^3.
    total = 0.0;
    for i = 1:n
        xi = i / n;
        for j = 1:n
            eta = j / n;
            for k = 1:n
                zeta = k / n;
                total = total + xi * eta^2 * zeta^2;
            end
        end
    end
end
