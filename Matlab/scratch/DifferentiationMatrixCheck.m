N = 10;
t0 = (1:10)';

D = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(1,t,1);
t1 = t0(2:end)-(t0(2)-t0(1))/2;

D12 = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(1,t1,1);
D1D2 = D12*D;

D2 = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(2,t,1);
t2 = t1(2:end)-(t1(2)-t1(1))/2;

D23 = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(1,t2,1);
D1D2D3 = D23*D1D2;

D3 = SmoothingSpline.FiniteDifferenceMatrixNoBoundary(3,t,1);