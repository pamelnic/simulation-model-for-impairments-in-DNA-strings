function [score, aligned_seq1, aligned_seq2, mutations] = needleman_wunsch( ...
        seq1, seq2, matchScore, mismatchScore, gapOpen, gapExtend, varargin)
% Gotoh NW with: asymmetric gaps, banded DP, auto-widen fallback, single-precision DP,
% O(n) traceback, and selectable logging ("full" | "counts").

    % Options 
    opts = struct;
    if ~isempty(varargin) && isstruct(varargin{1}), opts = varargin{1}; end
    if ~isfield(opts,'buildAlignmentStrings'), opts.buildAlignmentStrings = true; end
    if ~isfield(opts,'logMode'),               opts.logMode = 'full';     end
    if ~isfield(opts,'autoWiden'),             opts.autoWiden = true;     end

    % Asymmetric gaps
    if ~isfield(opts,'gapDelOpen'),   opts.gapDelOpen   = gapOpen;   end
    if ~isfield(opts,'gapDelExtend'), opts.gapDelExtend = gapExtend; end
    if ~isfield(opts,'gapInsOpen'),   opts.gapInsOpen   = gapOpen - 6; end
    if ~isfield(opts,'gapInsExtend'), opts.gapInsExtend = -2;          end

    % Inputs
    if isstring(seq1), seq1 = char(seq1); end
    if isstring(seq2), seq2 = char(seq2); end
    seq1 = upper(seq1(:).');  seq2 = upper(seq2(:).');
    n = length(seq1); m = length(seq2);

    % Band selection
    if ~isfield(opts,'bandW'), opts.bandW = []; end
    if ~isfield(opts,'maxBandW') || isempty(opts.maxBandW)
        opts.maxBandW = max(n,m) + 64;
    end

    % Cast scalars to single
    matchScore    = single(matchScore);
    mismatchScore = single(mismatchScore);
    gapDelOpen    = single(opts.gapDelOpen);
    gapDelExtend  = single(opts.gapDelExtend);
    gapInsOpen    = single(opts.gapInsOpen);
    gapInsExtend  = single(opts.gapInsExtend);

    % Internal worker so we can retry with wider band
    function [score_, a1_, a2_, mut_, endState] = solve_with_band(bw)
        N = n + 1; M = m + 1;

        H = ones(N, M, 'single') * (-inf('single'));
        E = ones(N, M, 'single') * (-inf('single'));
        F = ones(N, M, 'single') * (-inf('single'));

        Htb = uint8(zeros(N, M));  % 1=diag,2=fromE,3=fromF
        Etb = uint8(zeros(N, M));  % 0=open,1=extend
        Ftb = uint8(zeros(N, M));

        H(1,1) = single(0);
        for i = 2:N
            E(i,1)   = gapDelOpen + single(i-1) * gapDelExtend;
            H(i,1)   = E(i,1);
            Htb(i,1) = uint8(2); Etb(i,1) = uint8(1);
        end
        for j = 2:M
            F(1,j)   = gapInsOpen + single(j-1) * gapInsExtend;
            H(1,j)   = F(1,j);
            Htb(1,j) = uint8(3); Ftb(1,j) = uint8(1);
        end

        for i = 2:N
            ai = seq1(i-1);

            if isempty(bw)
                jStart = 2; jEnd = M;
            else
                jStart = max(2, i - bw);
                jEnd   = min(M, i + bw);
                if jStart > jEnd, continue; end
            end

            for j = jStart:jEnd
                bj = seq2(j-1);

                openE = H(i-1,j) + gapDelOpen   + gapDelExtend;
                extE  = E(i-1,j) + gapDelExtend;
                if openE >= extE
                    E(i,j) = openE; Etb(i,j) = uint8(0);
                else
                    E(i,j) = extE;  Etb(i,j) = uint8(1);
                end

                openF = H(i,  j-1) + gapInsOpen + gapInsExtend;
                extF  = F(i,  j-1) + gapInsExtend;
                if openF >= extF
                    F(i,j) = openF; Ftb(i,j) = uint8(0);
                else
                    F(i,j) = extF;  Ftb(i,j) = uint8(1);
                end

                s = matchScore; if ai ~= bj, s = mismatchScore; end
                diagScore = H(i-1,j-1) + s;

                if diagScore >= E(i,j) && diagScore >= F(i,j)
                    H(i,j) = diagScore; Htb(i,j) = uint8(1);
                elseif E(i,j) >= F(i,j)
                    H(i,j) = E(i,j);    Htb(i,j) = uint8(2);
                else
                    H(i,j) = F(i,j);    Htb(i,j) = uint8(3);
                end
            end
        end

        score_ = double(H(N,M));
        endState = Htb(N,M);

        % Traceback
        doStrings = opts.buildAlignmentStrings;

        % normalize logMode again 
        isFullLogs = (isstring(opts.logMode) && opts.logMode == "full") || ...
                     (ischar(opts.logMode)   && strcmpi(opts.logMode,'full'));
        wantFull  = isFullLogs;

        if doStrings
            maxLen = n + m;
            buf1 = repmat(' ', 1, maxLen);
            buf2 = repmat(' ', 1, maxLen);
            p = maxLen;
        else
            buf1 = []; buf2 = []; p = 0; 
        end

        if wantFull
            sub_log = struct('index', {}, 'original', {}, 'mutated', {});
            del_log = struct('index', {}, 'base_deleted', {});
            ins_log = struct('index_after', {}, 'base_inserted', {});
        else
            sub_log = struct('index', {}, 'original', {}, 'mutated', {});
            del_log = struct('index', {}, 'base_deleted', {});
            ins_log = struct('index_after', {}, 'base_inserted', {});
        end
        c_sub=0; c_del=0; c_ins=0;

        i = N; j = M;
        state = endState;
        if state == 0
            % band miss; caller may widen and retry
            a1_ = ''; a2_ = '';
            mut_.counts = struct('substitution',0,'deletion',0,'insertion',0);
            mut_.logs   = struct('substitution',sub_log,'deletion',del_log,'insertion',ins_log);
            return;
        end

        while (i > 1) || (j > 1)
            if i == 1 && j > 1, state = uint8(3); end
            if j == 1 && i > 1, state = uint8(2); end

            if state == 1
                a = seq1(i-1); b = seq2(j-1);
                if a ~= b
                    if wantFull
                        sub_log(end+1) = struct('index', i-1, 'original', a, 'mutated', b); 
                    end
                    c_sub = c_sub + 1;
                end
                if doStrings, buf1(p)=a; buf2(p)=b; p=p-1; end
                i = i - 1; j = j - 1; state = Htb(i,j);

            elseif state == 2
                a = seq1(i-1);
                if wantFull
                    del_log(end+1) = struct('index', i-1, 'base_deleted', a); 
                end
                c_del = c_del + 1;
                if doStrings, buf1(p)=a; buf2(p)='-'; p=p-1; end
                wasExt = (Etb(i,j) == uint8(1));
                i = i - 1;
                if wasExt
                    state = uint8(2);
                else
                    state = Htb(i,j);
                end

            else
                b = seq2(j-1);
                if wantFull
                    ins_log(end+1) = struct('index_after', i-1, 'base_inserted', b); 
                end
                c_ins = c_ins + 1;
                if doStrings, buf1(p)='-'; buf2(p)=b; p=p-1; end
                wasExt = (Ftb(i,j) == uint8(1));
                j = j - 1;
                if wasExt
                    state = uint8(3);
                else
                    state = Htb(i,j);
                end
            end
        end

        if doStrings
            a1_ = buf1(p+1:end); a2_ = buf2(p+1:end);
        else
            a1_ = ''; a2_ = '';
        end

        % Order left->right and compress insertions per line
        if wantFull
            if ~isempty(sub_log), sub_log = sub_log(end:-1:1); end
            if ~isempty(del_log), del_log = del_log(end:-1:1); end
            if ~isempty(ins_log), ins_log = ins_log(end:-1:1); end
            ins_log = merge_insertions(ins_log);
        end

        mut_.counts.substitution = c_sub;
        mut_.counts.deletion     = c_del;
        mut_.counts.insertion    = c_ins;
        mut_.logs.substitution   = sub_log;
        mut_.logs.deletion       = del_log;
        mut_.logs.insertion      = ins_log;

        aligned1 = a1_; aligned2 = a2_; mut = mut_; 
        a1_ = aligned1; a2_ = aligned2; mut_ = mut;
    end

    % Solve with band; widen on demand 
    curBand = opts.bandW;
    while true
        [score, aligned_seq1, aligned_seq2, mutations, endState] = solve_with_band(curBand);
        if endState ~= 0 || isempty(curBand) || ~opts.autoWiden
            break; % success, or full matrix, or no auto-widen
        end
        if curBand >= opts.maxBandW
            curBand = [];            % final attempt with full matrix
        else
            curBand = min(opts.maxBandW, max(curBand*2, curBand+16));
        end
    end
end