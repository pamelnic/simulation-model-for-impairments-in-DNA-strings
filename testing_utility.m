% Open files for reading line by line
originalFile = fopen("C:/Users/PC/MATLAB/Projects/Bachelor thesis_ver.1,3/sequences/dataset_testo1.txt", 'r');
mutatedFile = fopen("C:/Users/PC/MATLAB/Projects/Bachelor thesis_ver.1,3/sequences/dataset_testm3.txt", 'r');

totalBases = 0;
numInsertions = 0;
numDeletions = 0;
numSubstitutions = 0;

substitutionLog = {};
insertionLog = {};
deletionLog = {};

lineIndex = 1;
while ~feof(originalFile) && ~feof(mutatedFile)

    originalLine = fgetl(originalFile);
    mutatedLine = fgetl(mutatedFile);

    originalLine = regexprep(originalLine, '[^ACGT]', '');
    mutatedLine = regexprep(mutatedLine, '[^ACGT]', '');

    % Skip empty lines
    if isempty(originalLine) || isempty(mutatedLine)
        continue;
    end

    totalBases = totalBases + length(originalLine);

    % Needleman-Wunsch alignment
    [~, alignment] = nwalign(originalLine, mutatedLine, 'ALPHABET', 'NT');

    % Parse the alignment
    originalAligned = alignment(1, :);
    mutatedAligned = alignment(3, :);

    contextWindow = 2; 
    k = 1;
    while k <= length(originalAligned)
        if originalAligned(k) == '-'  % Insertion according to the logic
            insertionStart = k;
            while k < length(originalAligned) && originalAligned(k+1) == '-'
                k = k + 1;
            end
            insertionEnd = k;
           
            if insertionStart > 1 && insertionEnd < length(originalAligned) && ...
               originalAligned(insertionStart-1) == mutatedAligned(insertionEnd+1)
                % Insertion is likely due to misalignment, skip it
                k = k + 1;
                continue;
            end
            
            % log insertion
            numInsertions = numInsertions + 1;
            insertionLog{end+1} = sprintf('Line %d: Insertion from pos %d to %d', ...
                                          lineIndex, insertionStart, insertionEnd);

        elseif mutatedAligned(k) == '-'  % deletion
            numDeletions = numDeletions + 1;
            deletionLog{end+1} = sprintf('Line %d: Deletion at pos %d', lineIndex, k);

        elseif originalAligned(k) ~= mutatedAligned(k)  % maybe substitution
            startIndex = max(1, k - contextWindow);
            endIndex = min(length(originalAligned), k + contextWindow);

            % nogap near substitution
            if all(originalAligned(startIndex:endIndex) ~= '-') && ...
               all(mutatedAligned(startIndex:endIndex) ~= '-')
                numSubstitutions = numSubstitutions + 1;
                substitutionLog{end+1} = sprintf('Line %d: Substitution at pos %d (%s -> %s)', ...
                                                 lineIndex, k, originalAligned(k), mutatedAligned(k));
            end
        end
        k = k + 1;
    end

    lineIndex = lineIndex + 1;
end

fclose(originalFile);
fclose(mutatedFile);

% save logs
writecell(substitutionLog, 'C:\Users\PC\MATLAB\Projects\Bachelor Degree\logs/substitution_log.txt');
writecell(insertionLog, 'C:\Users\PC\MATLAB\Projects\Bachelor Degree\logs/insertion_log.txt');
writecell(deletionLog, 'C:\Users\PC\MATLAB\Projects\Bachelor Degree\logs/deletion_log.txt');

% Calculate rates
disp("Insertion count:-" + numInsertions);
disp("Subs:-" + numDeletions);
disp("Deletions:-" + numSubstitutions);
overall = ((numInsertions + numDeletions + numSubstitutions) / totalBases);
insertion = (numInsertions / totalBases);
deletion = (numDeletions / totalBases);
substitution = (numSubstitutions / totalBases);
errorRate = ((numInsertions + numDeletions + numSubstitutions) / totalBases) * 100;
insertionRate = (numInsertions / totalBases) * 100;
deletionRate = (numDeletions / totalBases) * 100;
substitutionRate = (numSubstitutions / totalBases) * 100;

fileSize = dir("C:/Users/PC/MATLAB/Projects/Bachelor thesis_ver.1,3/sequences/seqs_apollo_cleaned.txt").bytes

spreadsheetname= 'simulated_error_rates.xlsx';

if isfile(spreadsheetname)
    currentData = readcell(spreadsheetname);
    rowIndex = size(currentData, 1) + 1;
else
    rowIndex = 2;
    headers = {'Insertion', 'Deletion', 'Substitution', 'Overall'};
    writecell(headers, spreadsheetname, 'Range', 'A1');
end
probabilityVector = round([insertion, deletion, substitution, overall], 4);
modified = num2cell(probabilityVector);

writecell(modified, spreadsheetname, 'Range', sprintf('A%d', rowIndex));

% table with results(mostly for debugging purposes)
finalTable = table( ...
    fileSize, ... 
    lineIndex, ... 
    totalBases, ...
    errorRate, ...
    insertionRate, ...
    deletionRate, ...
    substitutionRate, ...
    'VariableNames', {'File_size_kB', 'Payloads', 'Bases', 'Error_rate', ...
                      'Insertion_rate', 'Deletion_rate', 'Substitution_rate'});
disp(finalTable);
disp('Logs saved to "logs/" directory.');
