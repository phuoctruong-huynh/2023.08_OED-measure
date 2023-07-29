% Compute d_HK(mu_1,mu_2)^2.
function distHK2 = compHK(q1, y1, q2, y2)

q1minus = (q1 < 0);
q1plus = (q1 > 0);
q2minus = (q2 < 0);
q2plus = (q2 > 0);

qq1 = [q1(q1plus); -q2(q2minus)];
q2 = [q2(q2plus); -q1(q1minus)];
q1 = qq1;
yy1 = [y1(q1plus); y2(q2minus)];
y2 = [y2(q2plus); y1(q1minus)];
y1 = yy1;

if isempty(y1) || isempty(y2)
  distHK2 = sum(q1) + sum(q2);
  return;
end

dist = abs(y1 - y2');

sindhs = min(sin(min(dist/2, pi/2)).^2, 1/2);
cosd = max(cos(min(dist, pi)), 0);

gamm = 1e-3 .* ones(size(cosd));

r1 = min(sqrt(q1 ./ sum(gamm, 2)), 1e30);
r2 = min(sqrt(q2' ./ sum(gamm, 1)), 1e30);

expgJh = cosd .* (r1 * r2);

for i = 1:1000

  %% multiplicative gradient update
  gamm = gamm .* expgJh;

  r1 = min(sqrt(q1 ./ sum(gamm, 2)), 1e30);
  r2 = min(sqrt(q2' ./ sum(gamm, 1)), 1e30);

  expgJh = cosd .* (r1 * r2);
end

ngamm = sum(gamm(:));

if ngamm > 0
  gamm = gamm / sum(gamm(:));
  gamm(gamm < 1e-15) = 0;
end

gamma1 = sum(gamm, 2);
gamma2 = sum(gamm, 1);

r1 = min(sqrt(q1 ./ gamma1), 1e30);
r2 = min(sqrt(q2' ./ gamma2), 1e30);

conic_d = (r1 - r2).^2 + 4 * (r1 * r2) .* sindhs;
distHK2 = sum(q1(gamma1==0)) + sum(q2(gamma2==0)) + gamm(:)' * conic_d(:);

end
