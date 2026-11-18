function result = f_generate_random_num_with_gap(range, count, decimal_places, min_gap)
% 生成指定范围内的随机数，任意两个数之间至少间隔 min_gap
%
% 输入参数：
%   range          - 范围，格式：[min, max]
%   count          - 生成的随机数个数
%   decimal_places - 保留的小数位数
%   min_gap        - 最小间隔（任意两个数的差值 >= min_gap）
%
% 输出：
%   result         - 生成的随机数向量（已排序）
%
% 示例：
%   generate_random_numbers_with_gap([0, 100], 5, 2, 10)
%   % 生成5个[0,100]之间的随机数，两两间隔至少10，保留2位小数

    % 参数检查
    if length(range) ~= 2 || range(1) >= range(2)
        error('范围必须为[min, max]格式，且min < max');
    end

    if ~isscalar(count) || count <= 0 || mod(count, 1) ~= 0
        error('count必须为正整数');
    end

    if ~isscalar(decimal_places) || decimal_places < 0 || mod(decimal_places, 1) ~= 0
        error('decimal_places必须为非负整数');
    end

    if ~isscalar(min_gap) || min_gap < 0
        error('min_gap必须为非负数');
    end

    % 检查是否可能生成满足条件的随机数
    min_val = range(1);
    max_val = range(2);
    total_range = max_val - min_val;

    min_span_needed = (count - 1) * min_gap;

    if min_span_needed > total_range
        error('范围[%.2f, %.2f]太小，无法放置%d个间隔至少%.2f的数', min_val, max_val, count, min_gap);
    end

    % 随机生成后检查，最多尝试10000次
    max_attempts = 10000;
    attempts = 0;

    while attempts < max_attempts
        % 生成随机数
        candidates = min_val + (max_val - min_val) * rand(count, 1);
        candidates = round(candidates, decimal_places);

        % 先排序再检查间隔；未排序直接 diff 会把随机顺序误当成空间顺序。
        candidates_sorted = sort(candidates);
        gaps = diff(candidates_sorted);

        if all(gaps >= min_gap)
            result_sort = candidates_sorted;
            indices = randperm(length(result_sort));
            result = result_sort(indices);
            return;
        end

        attempts = attempts + 1;
    end

    error('经过%d次尝试仍无法生成满足间隔要求的随机数，请增大范围或减小间隔', max_attempts);

end
