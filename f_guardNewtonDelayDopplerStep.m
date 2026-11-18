function guarded_step = f_guardNewtonDelayDopplerStep(raw_step, guard_threshold, max_guarded_step)
% 同时保护时延和多普勒 Newton 更新；两个 refinement 阶段共用本函数。
% 正常步长保持不变。任一维超过对应粗网格时，统一缩放二维步长，
% 保持 Newton 更新方向以及时延/多普勒之间的耦合比例。

    if any(~isfinite(raw_step))
        error('f_paraEstibyGradDesc_newtonOMP_multiTarget_v5:NonFiniteNewtonStep', ...
            'Newton delay/Doppler update contains NaN or Inf.');
    end

    delay_step_exceeded = abs(raw_step(1)) > guard_threshold(1);
    Doppler_step_exceeded = abs(raw_step(2)) > guard_threshold(2);
    guarded_step = raw_step;

    if delay_step_exceeded || Doppler_step_exceeded
        step_scale = max(abs(raw_step) ./ max_guarded_step);
        guarded_step = raw_step / step_scale;
        fprintf(['Newton异常步长保护已触发（时延超限=%d，多普勒超限=%d），' ...
            '缩放系数= %.4f\n'], delay_step_exceeded, Doppler_step_exceeded, step_scale);
    end
end