

write.csv(xcmsLargeWin@peaks, file = 'ori_peaks.csv')
write.csv(xcmsLargeWin@groups, file = 'ori_groups.csv')
ori_groupidx <- matrix(0,length(xcmsLargeWin@groupidx), max(unlist(lapply(xcmsLargeWin@groupidx, length))))
for (n in 1:length(xcmsLargeWin@groupidx))
{
  ori_groupidx[n,1:length(xcmsLargeWin@groupidx[[n]])] <- xcmsLargeWin@groupidx[[n]]
}
write.csv(ori_groupidx, file = 'ori_groupidx.csv')


write.csv(xcmsSmallWin@peaks, file = 'new_peaks.csv')
write.csv(xcmsSmallWin@groups, file = 'new_groups.csv')
new_groupidx <- matrix(0, length(xcmsSmallWin@groupidx), max(unlist(lapply(xcmsSmallWin@groupidx, length))))
for (n in 1:length(xcmsSmallWin@groupidx))
{
  new_groupidx[n, 1:length(xcmsSmallWin@groupidx[[n]])] <- xcmsSmallWin@groupidx[[n]]
}
write.csv(new_groupidx, file = 'new_groupidx.csv')





rt_raw <- matrix(0, max(sapply(agds@rt$raw, length)), length(agds@rt$raw))
for (n in 1:length(agds@rt$raw))
{
  temp_len <- length(agds@rt$raw[[n]])
  rt_raw[1:temp_len,n] <- agds@rt$raw[[n]]
}
write.csv(rt_raw, file = 'rt_raw.csv')


rt_cor2 <- matrix(0, max(sapply(agds@rt$corrected, length)), length(agds@rt$corrected))
for (n in 1:length(agds@rt$corrected))
{
  temp_len <- length(agds@rt$corrected[[n]])
  rt_cor2[1:temp_len,n] <- agds@rt$corrected[[n]]
}
write.csv(rt_cor2, file = 'rt_cor2.csv')


write.csv(data, file = 'data.csv')
