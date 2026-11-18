function y_query = sincInterpComplex(delay_axis, y, query_delay, half_len)
% sincInterpComplex
% 对复数距离匹配滤波结果进行窗化 sinc 插值。
%
% 输入：
%   delay_axis  - 匹配滤波结果对应的均匀 delay 网格，N_delay × 1
%   y           - 复数匹配滤波距离剖面，N_delay × 1
%   query_delay - 要读取的分数 delay，可以是任意列向量
%   half_len    - sinc 半长度，建议先取 8、12 或 16
%
% 输出：
%   y_query     - query_delay 处的复数插值结果
%
% 说明：
%   1. 直接插值复数结果，不能先取 abs 或 abs(.)^2；
%   2. 超出 delay_axis 的查询位置补零，不循环回绕；
%   3. 使用 Lanczos 窗截断 sinc。

    delay_axis = delay_axis(:);
    y = y(:);
    query_delay = query_delay(:);

    num_delay = numel(delay_axis);

    if num_delay < 2
        error('delay_axis 至少需要两个点。');
    end

    if numel(y) ~= num_delay
        error('y 的长度必须与 delay_axis 一致。');
    end

    delay_step = mean(diff(delay_axis));

    uniform_error = max(abs(diff(delay_axis) - delay_step));
    if uniform_error > 1e-9 * max(1, abs(delay_step))
        error('delay_axis 必须是均匀网格。');
    end

    % query_delay 在网格中的连续索引：
    % delay_axis(1) 对应 continuous_index = 1。
    continuous_index = ...
        (query_delay - delay_axis(1)) / delay_step + 1;

    center_index = round(continuous_index);
    tap_offset = -half_len : half_len;

    % 每一行对应一个查询位置。
    sample_index = center_index + tap_offset;

    valid_index = sample_index >= 1 & sample_index <= num_delay;
    inside_axis = continuous_index >= 1 & continuous_index <= num_delay;

    % 防止数组索引越界；无效位置的权重随后会被置零。
    safe_index = min(max(sample_index, 1), num_delay);

    % 查询点与各邻近采样点之间的归一化距离。
    distance = continuous_index - sample_index;

    % Lanczos 窗化 sinc。
    lanczos_order = half_len + 1;
    weight = sinc(distance) .* sinc(distance / lanczos_order);

    weight(~valid_index) = 0;

    % 归一化以减小有限长度截断产生的幅度误差。
    weight_sum = sum(weight, 2);
    good_query = inside_axis & abs(weight_sum) > 1e-12;

    weight(good_query, :) = ...
        weight(good_query, :) ./ weight_sum(good_query);

    weight(~good_query, :) = 0;

    % y(safe_index) 会形成“查询点 × sinc抽头”的复数矩阵。
    y_query = sum(weight .* y(safe_index), 2);
end