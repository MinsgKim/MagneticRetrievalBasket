% 파일에서 로그 데이터 읽기
filename = 'results.txt';
fileID = fopen(filename, 'r');
log_text = textscan(fileID, '%s', 'Delimiter', '\n');
fclose(fileID);

% 문자열 배열로 변환
log_lines = string(log_text{1});

% 정규표현식으로 Bx, By, Bz, Norm 추출 (음수, 지수도 포함)
pattern = 'Bx = ([\-0-9.eE]+), By = ([\-0-9.eE]+), Bz = ([\-0-9.eE]+), Norm = ([\-0-9.eE]+)';
tokens = regexp(log_lines, pattern, 'tokens');

% 유효한 라인만 필터링
valid_idx = ~cellfun(@isempty, tokens);
tokens = tokens(valid_idx);

% 추출된 문자열을 숫자로 변환
numPoints = length(tokens);
Bx = zeros(numPoints,1);
By = zeros(numPoints,1);
Bz = zeros(numPoints,1);
Norm = zeros(numPoints,1);

for i = 1:numPoints
    temp = tokens{i}{1};  % 내부 셀 배열 가져오기
    Bx(i) = str2double(temp{1});
    By(i) = str2double(temp{2});
    Bz(i) = str2double(temp{3});
    Norm(i) = str2double(temp{4});
end

% 테이블로 정리
T = table(Bx, By, Bz, Norm);

% Excel로 저장
writetable(T, 'cubic_magnetic_field_data2.xlsx');

disp('Excel 파일 저장 완료: magnetic_field_data.xlsx');
