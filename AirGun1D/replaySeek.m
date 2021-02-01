tSeek = .1545;
tCutoff = .1652;

indexFirst = find(solution.soln.x > tSeek, 1, 'first');
indexLast = find(solution.soln.x < tCutoff, 1, 'last');
