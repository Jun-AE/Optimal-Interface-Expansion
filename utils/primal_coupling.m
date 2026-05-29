function [L, K_AB, M_AB, C_AB] = primal_coupling(ninterface, KA, KB, MA, MB, CA, CB)
% PRIMAL_COUPLING  Primally couple two substructures across a shared interface.
%
%   [L, K_AB, M_AB, C_AB] = PRIMAL_COUPLING(NINTERFACE, KA, KB, MA, MB, CA, CB)
%   couples substructures A and B by enforcing compatibility at NINTERFACE
%   shared interface nodes. The compatible (reduced) DoF vector q relates to
%   the uncoupled block-diagonal vector u via u = L*q, and the coupled
%   system matrices follow:
%
%       K_AB = L' * blkdiag(KA, KB) * L
%       M_AB = L' * blkdiag(MA, MB) * L
%       C_AB = L' * blkdiag(CA, CB) * L
%
%   This is the explicit primal projection from LM-FBS (paper §2.2, Eqs. 6-7).
%   L is built directly rather than via NULL(B); the construction below is
%   equivalent to L = null(B) for the layout described in the contract.
%
%   Inputs:
%     NINTERFACE - (integer, scalar) Number of shared interface NODES (not
%                  DoFs). For 2-DoF/node Euler-Bernoulli beam elements the
%                  number of interface DoFs is 2*NINTERFACE.
%     KA, KB     - (nA×nA, nB×nB double) Substructure stiffness matrices.
%     MA, MB     - (nA×nA, nB×nB double) Substructure mass matrices.
%     CA, CB     - (nA×nA, nB×nB double) Substructure damping matrices.
%
%   Outputs:
%     L    - ((nA+nB) × (nA+nB-2*NINTERFACE) double) Primal projector u = L*q.
%            Currently unused by the calling script but returned for API
%            stability and downstream diagnostic use.
%     K_AB - ((nA+nB-2*NINTERFACE) × (nA+nB-2*NINTERFACE)) Coupled stiffness.
%     M_AB - same size as K_AB. Coupled mass.
%     C_AB - same size as K_AB. Coupled damping.
%
%   DoF layout contract (precondition — silently violated if not met):
%     - Substructure A's interface DoFs are its LAST  2*NINTERFACE rows/cols.
%     - Substructure B's interface DoFs are its FIRST 2*NINTERFACE rows/cols.
%     - Interface DoFs are paired by STRAIGHT one-to-one correspondence:
%
%         A.dof(end - 2*ninterface + k + 1)  ==  B.dof(k + 1)
%         for k = 0 .. 2*ninterface - 1
%
%       i.e. A's first interface DoF (lowest index inside its interface
%       band) is identified with B's first interface DoF, A's second with
%       B's second, and so on. This requires the caller to number A's and
%       B's interface DoFs in the same ordinal order — the natural case
%       when both substructures number their nodes from one end to the
%       other along the coupling direction.
%
%   Element-type assumption:
%     The 2*NINTERFACE conversion assumes 2 DoFs per node (transverse
%     displacement + rotation), as in the project's Euler-Bernoulli beam
%     model. For shell or 3D elements, pass the equivalent node count or
%     adapt this function.
%
%   Example:
%     % Couple two cantilever beams sharing 3 interface nodes (= 6 DoFs):
%     [L, Kfull, Mfull, Cfull] = primal_coupling(3, ...
%         beamA_exp.K, beamB_exp.K, ...
%         beamA_exp.M, beamB_exp.M, ...
%         beamA_exp.C, beamB_exp.C);
%
%   Reference:
%     Junaid et al. (2026). Journal of Sound and Vibration. §2.2 (LM-FBS).
%     DOI: 10.1016/j.jsv.2026.001458
%     de Klerk, D., Rixen, D.J., & Voormeeren, S.N. (2008). General framework
%       for dynamic substructuring. AIAA Journal 46(5): 1169-1181.
%
%   See also: couple_substructures, blkdiag, eig.

%% Interface DoF count

% Convert interface-node count to interface-DoF count.
% Assumes 2 DoFs/node (Euler-Bernoulli beam); change here for other elements.
nInterfaceDoFs = 2 * ninterface;

validate_interface_dofs(ninterface, size(KA, 1), size(KB, 1), ...
    'Layout', 'full', 'Caller', 'primal_coupling');

% Substructure sizes and total uncoupled DoF count.
nA      = size(KA, 1);
nB      = size(KB, 1);
ntdofs  = nA + nB;                       % total DoFs in uncoupled blocks
nq      = ntdofs - nInterfaceDoFs;       % DoFs in the coupled (reduced) set

%% L matrix construction

% L1 maps q -> uA. Identity in the first nA rows; the trailing
% nInterfaceDoFs columns of L1 represent the shared interface band.
L1 = eye(nA, nq);

% L2 maps q -> uB. rot90(L1, 2) (i.e. flipud+fliplr) shifts the identity
% block so that L2's leading nInterfaceDoFs rows pick up the SAME q
% indices that L1's trailing rows did — yielding the straight 1-to-1
% interface pairing described in the docstring contract.
L2 = rot90(L1, 2);

L = [L1; L2];

%% Block-diagonal global matrices (uncoupled system)

K_combined = blkdiag(KA, KB);
M_combined = blkdiag(MA, MB);
C_combined = blkdiag(CA, CB);

%% Primal projection: u = L*q  =>  L'*K*L

K_AB = L' * K_combined * L;
M_AB = L' * M_combined * L;
C_AB = L' * C_combined * L;

end
