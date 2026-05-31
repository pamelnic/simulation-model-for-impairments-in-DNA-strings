
function [summary_del, summary_ins, summary_sub] = iid_channel% Open files for reading and writing sequences
    fileRID = fopen('/home/paulpc/Projects/Bachelor thesis_ver 13/Bachelor thesis_ver.1,3/sequences/seqs_apollo_cleaned.txt','r');
    fileWID = fopen('/home/paulpc/Projects/Bachelor thesis_ver 13/Bachelor thesis_ver.1,3/sequences/mut-seq_apollo.txt','w');
    if fileRID == -1 || fileWID == -1
        if fileRID == -1
            warning('Failed to open input file: sequences\\seqs_apollo_cleaned.txt');
        end
        if fileWID == -1
            warning('Failed to open output file: sequences\\mut-seq_apollo.txt');
        end
        if fileRID == -1 && fileWID ~= -1
            fclose(fileWID);
        end
        error('Cannot continue: input or output file could not be opened.');
    end
    
    % Probabilities for deletion, insertion, and substitution
    p_del = 0.0295;
    p_ins = 0.0216;
    p_substitution = 0.0175;
    
    % Substitution probability matrix and nucleotides
    sub_matrix = [0,     0.043, 0.190, 0.033;  
                  0.056, 0,     0.025, 0.130;  
                  0.120, 0.023, 0,     0.047;  
                  0.035, 0.260, 0.030, 0    ]; 
    nucleotides = 'ACGT';
    
    % Counters testing
    n_del = 0; n_ins = 0; n_sub = 0;
    orig_bases = 0; mut_bases_total = 0;
    sub_counts = zeros(4,4); % how many times i->j happened

    


    line = fgetl(fileRID);
    while ischar(line)
        % clean a bit: uppercase
        line = upper(line);
        orig_bases = orig_bases + numel(line);
    
        mut_sequence = line;
        i = 1;
        while i <= length(mut_sequence)
            r = rand();
    
            if r < p_del
                % Deletion: remove current char; do not advance i
                mut_sequence(i) = [];
                n_del = n_del + 1;
    
            elseif r < p_del + p_ins
                % Insertion: insert BEFORE current char; skip the inserted char
                new_nucleotide = nucleotides(randi(4));
                if i == 1
                    mut_sequence = [new_nucleotide, mut_sequence];
                else
                    mut_sequence = [mut_sequence(1:i-1), new_nucleotide, mut_sequence(i:end)];
                end
                n_ins = n_ins + 1;
                i = i + 1; % skip the inserted char, process the original current next
    
            elseif r < p_del + p_ins + p_substitution
                % Substitution: only if current is ATCG
                current_nucleotide = mut_sequence(i);
                idx = find(nucleotides == current_nucleotide, 1);
                if ~isempty(idx)
                    weights = sub_matrix(idx, :);
                    % Draw with weights; avoid randsample dependency if missing:
                    new_nucleotide = weighted_pick(nucleotides, weights);
                    mut_sequence(i) = new_nucleotide;
    
                    % log substitution matrix counts
                    jdx = find(nucleotides == new_nucleotide, 1);
                    if ~isempty(jdx), sub_counts(idx, jdx) = sub_counts(idx, jdx) + 1; end
    
                    n_sub = n_sub + 1;
                end
                i = i + 1;
    
            else
                i = i + 1;
            end
        end
    
        fprintf(fileWID, '%s\n', mut_sequence);
        mut_bases_total = mut_bases_total + numel(mut_sequence);
    
        line = fgetl(fileRID);
    end
    
    % Close files
    fclose(fileRID);
    fclose(fileWID);
    
    % Summary
    summary_del = n_del/max(1,orig_bases);
    summary_sub = n_sub/max(1,orig_bases);
    summary_ins = n_ins/max(1,orig_bases);
   
    fprintf('\n--- Mutation summary ---\n');
    fprintf('Original bases (sum of line lengths): %d\n', orig_bases);
    fprintf('Mutated  bases (sum of line lengths): %d\n', mut_bases_total);
    fprintf('Deletions:    %d  (%.4f per original base)\n', n_del, n_del/max(1,orig_bases));
    fprintf('Insertions:   %d  (%.4f per original base)\n', n_ins, n_ins/max(1,orig_bases));
    fprintf('Substitutions:%d  (%.4f per original base)\n', n_sub, n_sub/max(1,orig_bases));

     % Display substitution count matrix
    row_sums_ref = sum(sub_matrix, 2);
    sub_matrix_norm = sub_matrix ./ max(row_sums_ref, eps);
    
    fprintf('\n--- Normalized reference substitution matrix ---\n');
    fprintf('Rows = original nucleotide, Columns = substituted nucleotide\n\n');
    
    disp(array2table(sub_matrix_norm, ...
        'VariableNames', cellstr(nucleotides'), ...
        'RowNames', cellstr(nucleotides')));

    % Display substitution probabilities estimated from simulation
    row_sums = sum(sub_counts, 2);
    sub_freq = sub_counts ./ max(row_sums, 1);

    fprintf('\n--- Simulated substitution frequency matrix ---\n');
    fprintf('Each row shows how often one nucleotide was replaced by another\n\n');

    disp(array2table(sub_freq, ...
        'VariableNames', cellstr(nucleotides'), ...
        'RowNames', cellstr(nucleotides')));
    
    if any(sum(sub_matrix,2)==0)
        warning('One or more substitution rows sum to zero; substitution drawing would fail.');
    end

    function ch = weighted_pick(pop, weights)
        w = double(weights(:)');
        s = sum(w);
        if s <= 0
            % fallback to uniform over nonzero weights; if all zero, uniform over all
            mask = w > 0;
            if any(mask), w = mask / sum(mask); else, w = ones(size(w)) / numel(w); end
        else
            w = w / s;
        end
        u = rand();
        k = find(u <= cumsum(w), 1, 'first');
        ch = pop(k);
    end

end