dat_P = monitorStates(13,:);
dat_t = monitorStates(14,:);

av_P = 0.5 * (dat_P(1:end-1) + dat_P(2:end));
av_t = 0.5 * (dat_t(1:end-1) + dat_t(2:end));
dat_dt = diff(dat_t);

delta_Energy = nan(length(av_P),1);
for i = 1:length(dat_dt)
    delta_Energy(i) = sum(av_P(1:i) .* dat_dt(1:i));
end
plot(delta_Energy);