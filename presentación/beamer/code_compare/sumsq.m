% =============================================================================
%  sumsq.m — Evaluación de una función de base hexaédrica de Nédélec
%  (componente xi de H2) sobre una nube de puntos del hexaedro de referencia.
%
%  MATLAB/Octave NO permiten usar caracteres griegos como identificadores,
%  por lo que hay que escribir xi, eta, zeta con nombres ASCII.
%
%  Ejecutar:  octave sumsq.m   (o  matlab -batch sumsq)
% =============================================================================

n = 200;
tic;
result = evaluate(n);
t = toc;
fprintf('MATLAB  n=%d  resultado=%g  tiempo=%.2f ms\n', n, result, t*1000);
