% DNA alignment script using parallel processing
clear; % Clear all variables

opts = struct('bandW',[], 'autoWiden',false, ...
              'buildAlignmentStrings',false, ...
              'logMode','counts','countMode','base', ...
              'gapInsOpen', -27, 'gapInsExtend', -3);

% Test: one deleted base at the start
[~,~,~,m1] = needleman_wunsch('GACGT','ACGT', 5,-19,-19,-1, opts);
disp(m1.counts)   % expected: 1 deletion

% Test: one inserted base at the end
[~,~,~,m2] = needleman_wunsch('ACGT','ACGTT', 5,-19,-19,-1, opts);
disp(m2.counts)   % expected: 1 insertion


% 1) Input files 
original_file = '/home/paulpc/Projects/Bachelor thesis_ver 13/Bachelor thesis_ver.1,3/sequences/seqs_apollo_cleaned.txt';
mutated_file  = '/home/paulpc/Projects/Bachelor thesis_ver 13/Bachelor thesis_ver.1,3/sequences/mut-seq_apollo.txt';

% 2) Alignment scoring 
matchScore     =  5;
mismatchScore  = -24;
gapOpen        = -23; 
gapExtend      =  -1;

% 3) Program settings 
DO_FULL_LOGS            = false;  % true = save detailed mutation logs
BUILD_ALIGNMENT_STRINGS = false;   % true = store aligned sequences

% 4) Needleman-Wunsch options 
nwOpts = struct();
nwOpts.bandW                 = 24;    % starting band width
nwOpts.autoWiden             = true;  % increase band if needed
nwOpts.maxBandW              = 256;   % maximum allowed band width
nwOpts.buildAlignmentStrings = BUILD_ALIGNMENT_STRINGS;

if DO_FULL_LOGS
    nwOpts.logMode = 'full';
else
    nwOpts.logMode = 'counts';
end

% Make insertions slightly harder to select
nwOpts.gapInsOpen            = gapOpen - 5;
nwOpts.gapInsExtend          = -3;

% 5) Start parallel pool 
pool = gcp('nocreate');

if isempty(pool)
    try
        parpool("Threads");   % use thread pool if available
    catch
        parpool;              % otherwise use normal process pool
    end
end

% 6) Run alignment
disp(' Starting Parallel Segmented Alignment ');
tic;

try
    [final_score, aligned_orig, aligned_mut, mutations] = run_parallel_alignment_fast( ...
        original_file, mutated_file, matchScore, mismatchScore, gapOpen, gapExtend, nwOpts);

    t = toc;

    % Total number of mutation events used for percentage calculation
    total_mutations = 15312600;

    if total_mutations > 0
        perc_sub = 100 * double(mutations.counts.substitution) / total_mutations;
        perc_del = 100 * double(mutations.counts.deletion)     / total_mutations;
        perc_ins = 100 * double(mutations.counts.insertion)    / total_mutations;
    else
        [perc_sub, perc_del, perc_ins] = deal(0);
    end

    % Display final results
    disp('---------------------------------');
    disp(['Processing Time:      ', num2str(t, '%.2f'), ' seconds']);
    disp(['Optimal Score (sum):  ', num2str(final_score)]);
    disp('---------------------------------');
    disp('--- Mutation Summary (Percentage Breakdown) ---');
    disp(['Total Mutation Events: ', num2str(total_mutations)]);
    disp(['Substitutions: ', num2str(mutations.counts.substitution), ' (', num2str(perc_sub, '%.2f'), '%)']);
    disp(['Deletions:     ', num2str(mutations.counts.deletion),     ' (', num2str(perc_del, '%.2f'), '%)']);
    disp(['Insertions:    ', num2str(mutations.counts.insertion),    ' (', num2str(perc_ins, '%.2f'), '%)']);
    disp('----------------------------------------------');

catch ME
    t = toc;

    % Show error message if alignment fails
    fprintf('Error during parallel alignment execution (after %.2f s):\n', t);
    fprintf('%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));

    delete(gcp('nocreate'));
end


%% LOCAL FUNCTIONS 

function [final_score, aligned_orig, aligned_mut, mutations] = run_parallel_alignment_fast( ...
    file1, file2, matchScore, mismatchScore, gapOpen, gapExtend, opts)

% Runs alignment for all sequence lines in parallel.


    disp('Step 1/3: Reading all lines into memory...');
    seq1_lines = read_all_lines(file1);
    seq2_lines = read_all_lines(file2);

    L1 = numel(seq1_lines);
    L2 = numel(seq2_lines);
    L  = max(L1, L2);

    disp(['Step 2/3: Aligning ', num2str(L), ' line pairs...']);

    % Preallocate arrays for speed
    scores  = zeros(L, 1);
    counts  = zeros(L, 3, 'uint32');  % columns: substitutions, deletions, insertions

    wantStrings = isfield(opts,'buildAlignmentStrings') && opts.buildAlignmentStrings;

    % Check if detailed logs are needed
    if ~isfield(opts,'logMode')
        opts.logMode = 'full';
    end

    isFullLogs = (isstring(opts.logMode) && opts.logMode == "full") || ...
                 (ischar(opts.logMode)   && strcmpi(opts.logMode,'full'));

    wantFull = isFullLogs;

    % Prepare storage for aligned strings if needed
    if wantStrings
        aligned_orig_segments = cell(L, 1);
        aligned_mut_segments  = cell(L, 1);
    else
        aligned_orig_segments = [];
        aligned_mut_segments  = [];
    end

    % Prepare storage for full logs if needed
    if wantFull
        seg_logs_sub = cell(L,1);
        seg_logs_del = cell(L,1);
        seg_logs_ins = cell(L,1);
    end

    % Align each pair of lines in parallel
    parfor line_number = 1:L
        s1 = '';

        if line_number <= L1
            s1 = seq1_lines{line_number};
        end

        s2 = '';

        if line_number <= L2
            s2 = seq2_lines{line_number};
        end

        [sc, a1, a2, segMut] = needleman_wunsch( ...
            s1, s2, matchScore, mismatchScore, gapOpen, gapExtend, opts);

        scores(line_number) = sc;

        % Save mutation counts for this line
        counts(line_number, :) = uint32([ ...
            segMut.counts.substitution, segMut.counts.deletion, segMut.counts.insertion]);

        if wantStrings
            aligned_orig_segments{line_number} = a1;
            aligned_mut_segments{line_number}  = a2;
        end

        if wantFull
            seg_logs_sub{line_number} = segMut.logs.substitution;
            seg_logs_del{line_number} = segMut.logs.deletion;
            seg_logs_ins{line_number} = segMut.logs.insertion;
        end
    end

    disp('Step 3/3: Combining results...');

    % Sum all alignment scores
    final_score = sum(scores);

    % Join aligned strings only if they were created
    if wantStrings
        aligned_orig = strjoin(aligned_orig_segments, newline);
        aligned_mut  = strjoin(aligned_mut_segments,  newline);
    else
        aligned_orig = '';
        aligned_mut  = '';
    end

    % Calculate total mutation counts
    mutations.counts.substitution = sum(counts(:,1), 'native');
    mutations.counts.deletion     = sum(counts(:,2), 'native');
    mutations.counts.insertion    = sum(counts(:,3), 'native');

    % Store per-line mutation counts
    mutations.per_line_counts = counts;

    % Combine detailed logs if full logging is enabled
    if wantFull
        for k = 1:L
            s = seg_logs_sub{k};

            if ~isempty(s)
                [s.line_number] = deal(k);
                seg_logs_sub{k} = s;
            end

            d = seg_logs_del{k};

            if ~isempty(d)
                [d.line_number] = deal(k);
                seg_logs_del{k} = d;
            end

            z = seg_logs_ins{k};

            if ~isempty(z)
                [z.line_number] = deal(k);
                seg_logs_ins{k} = z;
            end
        end

        subs = vertcat_safe(seg_logs_sub);
        dels = vertcat_safe(seg_logs_del);
        ins  = vertcat_safe(seg_logs_ins);

        % Merge nearby insertions into cleaner events
        ins = merge_insertions(ins);

        mutations.logs.substitution = subs;
        mutations.logs.deletion     = dels;
        mutations.logs.insertion    = ins;

        % Update counts using the final logs
        mutations.counts.substitution = uint32(numel(subs));
        mutations.counts.deletion     = uint32(numel(dels));
        mutations.counts.insertion    = uint32(numel(ins));
    else
        % Empty logs when only counts are used
        mutations.logs.substitution = struct('line_number',{},'index',{},'original',{},'mutated',{});
        mutations.logs.deletion     = struct('line_number',{},'index',{},'base_deleted',{});
        mutations.logs.insertion    = struct('line_number',{},'index_after',{},'base_inserted',{});
    end
end


function out = vertcat_safe(cells)
% Combines non-empty struct arrays from a cell array.

    if isempty(cells)
        out = struct([]);
        return;
    end

    nonempty = ~cellfun('isempty', cells);

    if ~any(nonempty)
        out = struct([]);
        return;
    end

    out = vertcat(cells{nonempty});
end


function line_array = read_all_lines(filename)
% Reads all lines from a text file and converts them to uppercase.

    fid = fopen(filename, 'r');

    if fid == -1
        error('Cannot open file: %s', filename);
    end

    C = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    fclose(fid);

    line_array = C{1};

    % Remove extra spaces and make all letters uppercase
    line_array = cellfun(@(x) upper(strtrim(x)), line_array, 'UniformOutput', false);
end