% ===========================================================================
% Benchmark de rendimiento ensamblaje y solver MATLAB
% ===========================================================================
%
% Etapas (comparables 1:1 con benchmark_julia.jl):
%   coefficients   getCoefficients()  
%   dof_numbering  mesh.generateDOF() 
%   assembly       assembly.getGlobalMatrices() — matrices globales K y M.
%                  en MATLAB esta función también impone las condiciones
%                  Dirichlet (imposeLocalDirichletMatrix) y en julia se hace después
%   solver         eigs(K, M, numberKc, sigma) — los 6 modos físicos más bajos,
%                  
%
% Metodología:
%   - 1 Hilo
%   - TTFX (primera llamada, tic/toc) y timeit que hace múltiples repeticiones
%
% Mallas: projectMeshes de meshesHalfFilledTAP2021.mat.
%
% Salida: output/results_matlab.csv .
% ===========================================================================

% --- Localizar rutas relativas a este script -------------------------------
thisFile   = mfilename('fullpath');
testDir    = fileparts(thisFile);
matlabDir  = fullfile(testDir, '..', '..', 'matlab');
outputDir  = fullfile(testDir, '..', 'output');
if ~exist(outputDir, 'dir'); mkdir(outputDir); end
csvPath    = fullfile(outputDir, 'results_matlab.csv');

% getCoefficients / getGlobalMatrices usan rutas relativas a matlab/ (+bases,
% +settings…), 
origDir = pwd;
cd(matlabDir);

% Un solo hilo de cómputo: comparación de algoritmo serial, no de paralelismo.
maxNumCompThreads(1);

% Las bases precomputadas se renombran a .bak 
% y se restauran al final para forzar el calcularlas sin borrar las antiguas
precomputedFile = fullfile('+bases', 'precomputedBases', 'systematicBasesHexa2.mat');
backupFile      = [precomputedFile '.bak'];
if isfile(precomputedFile); movefile(precomputedFile, backupFile); end

fid = -1;
try
    % --- Fichero de salida con cabecera común ------------------------------
    fid = fopen(csvPath, 'w');
    fprintf(fid, ['lenguaje,malla,n_elem,n_dof,etapa,t_ttfx_ms,t_warm_ms,' ...
                  'bytes_alloc,mem_occupied_bytes,nnz,nnz_free,mem_free_bytes\n']);

    % -------------------------------------------------------------------
    % coeficientes 
    % -------------------------------------------------------------------
    disp('[matlab] etapa coefficients ...');
    sysBasesObj = bases.SystematicBases();
    tic; sysBasesObj.getCoefficients(); tTtfx = toc * 1e3;     % primera vez: real
    tWarm = timeit(@() sysBasesObj.getCoefficients()) * 1e3;   % repeticiones
    emitRow(fid, 0, 0, 0, 'coefficients', tTtfx, tWarm, NaN, NaN, NaN);
    sysBasesObj.saveCoefficients();   % dejarlas en la caché para el ensamblaje

    % -------------------------------------------------------------------
    % dof_numbering, assembly, solver
    % -------------------------------------------------------------------
    load(fullfile('+settings', 'meshesHalfFilledTAP2021.mat'), 'projectMeshes');
    nMeshes = length(projectMeshes);

    % Autovalores analíticos (rad/m) de la cavidad semillena, para centrar sigma de
    % eigs en la banda física, igual que mainHalfFilled.m. 
    kcAnalytic = [353.78 544.50 599.80 750.31 763.32 809.58];
    numberKc   = 6;
    sigmaEst   = kcAnalytic(1)^2 + 1.5 * (kcAnalytic(end) - kcAnalytic(1))^2;

    for i = 1:nMeshes
        meshObj = projectMeshes(i);
        nElem   = length(meshObj.elements);
        fprintf('[matlab] malla %d (%d elem) ...\n', i, nElem);

        % --- dof_numbering -------------------------------------------------
        % generateDOF muta el objeto (Mesh es un handle class),llamarlo varias veces hace 
        % lo mismo que llamarlo una por lo que si lo llamas varias veces hay que reiniciar currentDOF =0 
        % para que se sobreescriba y no se acumule.
        tic; nDof = meshObj.generateDOF(); tTtfx = toc * 1e3;
        tWarm = timeit(@() meshObj.generateDOF()) * 1e3;
        emitRow(fid, i, nElem, nDof, 'dof_numbering', tTtfx, tWarm, NaN, NaN, NaN);

        % --- assembly ------------------------------------------------------
        % evalc suprime el "X% assembly done" para no ensuciar la consola ni el cronómetro .
        tic; [~] = evalc('FEM = assembly.getGlobalMatrices(meshObj);'); tTtfx = toc * 1e3;
        tWarm = timeit(@() suppressAssembly(meshObj)) * 1e3;
        K = FEM.stiffness; M = FEM.mass;
        sK = whos('K'); sM = whos('M');
        memOcc = sK.bytes + sM.bytes;

        % Bloque libre-libre. La K de MATLAB es n_dof x n_dof y la de Julia es n_free x n_free, 
        %Lo comparable en ambos lados es K(free,free).
        freeDOF = getFreeDOF(meshObj, nDof);
        Kf = K(freeDOF, freeDOF); Mf = M(freeDOF, freeDOF);
        sKf = whos('Kf'); sMf = whos('Mf');
        memFree = sKf.bytes + sMf.bytes;
        emitRow(fid, i, nElem, nDof, 'assembly', tTtfx, tWarm, NaN, memOcc, nnz(K), ...
                nnz(Kf), memFree);

        % --- solver (autovalores) ------------------------------------------
        % Los 6 modos físicos más bajos con eigs shift-invert como mainHalfFilled.m:63-64.
        tic; [~] = evalc('eigs(K, M, numberKc, sigmaEst);'); tTtfx = toc * 1e3;
        tWarm = timeit(@() eigsQuiet(K, M, numberKc, sigmaEst)) * 1e3;
        emitRow(fid, i, nElem, nDof, 'solver', tTtfx, tWarm, NaN, NaN, NaN);
    end

    cleanupBench(fid, origDir, precomputedFile, backupFile);
    fprintf('[matlab] listo -> %s\n', csvPath);
catch ME
    % Cerrar el CSV, restaurar la caché de bases y volver al directorio original
    % antes de propagar el error: un fallo a mitad no debe dejar residuos.
    cleanupBench(fid, origDir, precomputedFile, backupFile);
    rethrow(ME);
end

% ===========================================================================
% Funciones auxiliares
% ===========================================================================

% Escribe una fila del CSV con el esquema común. Los campos NaN se imprimen "NA".
% nnzFree/memFree (opcionales) se refieren al bloque libre-libre.
function emitRow(fid, malla, n_elem, n_dof, etapa, ttfxMs, warmMs, ...
                 bytesAlloc, memOcc, nnzVal, nnzFree, memFree)
    if nargin < 11; nnzFree = NaN; end
    if nargin < 12; memFree = NaN; end
    fprintf(fid, 'matlab,%d,%d,%d,%s,%.6f,%.6f,%s,%s,%s,%s,%s\n', ...
            malla, n_elem, n_dof, etapa, ttfxMs, warmMs, ...
            fld(bytesAlloc), fld(memOcc), fld(nnzVal), fld(nnzFree), fld(memFree));
end

% DOFs libres = todos menos los de aristas/caras en una polyFace "PerfectE".
% Es la MISMA condición que aplica imposeLocalDirichletMatrix dentro del
% ensamblaje (y la misma que usa exportGlobalMatrices.m).
function freeDOF = getFreeDOF(meshObj, nDof)
    dirichletDOF = [];
    for ie = 1:length(meshObj.edges)
        pf = meshObj.edges(ie).polyFaceId;
        if (pf ~= 0 && strcmp(meshObj.polyFaces(pf).polyFaceType, "PerfectE"))
            dirichletDOF = [dirichletDOF; meshObj.edges(ie).dof(:)]; %#ok<AGROW>
        end
    end
    for iface = 1:length(meshObj.faces)
        pf = meshObj.faces(iface).polyFaceId;
        if (pf ~= 0 && strcmp(meshObj.polyFaces(pf).polyFaceType, "PerfectE"))
            dirichletDOF = [dirichletDOF; meshObj.faces(iface).dof(:)]; %#ok<AGROW>
        end
    end
    freeDOF = setdiff((1:nDof)', unique(dirichletDOF));
end

% Formatea un campo numérico; NaN -> "NA".
function s = fld(x)
    if isnan(x); s = 'NA'; else; s = num2str(x); end
end

% Ejecuta el ensamblaje suprimiendo su salida por pantalla, para timeit.
function suppressAssembly(meshObj)
    [~] = evalc('assembly.getGlobalMatrices(meshObj);');
end

% Ejecuta eigs suprimiendo posibles avisos, para timeit.
function eigsQuiet(K, M, numberKc, sigmaEst)
    [~] = evalc('eigs(K, M, numberKc, sigmaEst);');
end

% Deja el entorno como estaba: CSV cerrado, caché de bases restaurada y
% directorio de trabajo original. Es idempotente, para poder llamarla tanto en
% la salida normal como en el catch.
function cleanupBench(fid, origDir, precomputedFile, backupFile)
    if fid > 0; fclose(fid); end
    if isfile(backupFile)
        if isfile(precomputedFile)
            delete(backupFile);            % se regeneró: el .bak sobra
        else
            movefile(backupFile, precomputedFile);   % no se regeneró: restaurar
        end
    end
    cd(origDir);
end
