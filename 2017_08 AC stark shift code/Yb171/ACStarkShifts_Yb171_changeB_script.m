% script to see how magnitude of magnetic affects AC shifts.
% conclusion: changing B from 20 to 50 has ~2% increase on magnitude of
% absolute Stark shift

B = 0 : 5 : 100;
diff_1 = zeros(size(B));
diff_end = diff_1;
shift_up = diff_1;
shift_down = diff_1;

for i = 1 : length(B)
    b = B(i);
    [AC_diff_1, AC_diff_end, abs_shift_up, abs_shift_down] = ACStarkShifts_Yb171_changeB( b );
    diff_1(i) = AC_diff_1;
    diff_end(i) = AC_diff_end;
    shift_up(i) = abs_shift_up;
    shift_down(i) = abs_shift_down;
end

figure(1), clf, hold on
plot(B, diff_1, 'b-*')
plot(B, diff_end, 'r-*')
hold off
xlabel('B /gauss'), ylabel('Diff AC Stark shift /Hz')
legend('Diff at start', 'Diff at end')
hold off

figure(2), clf, hold on
plot(B, shift_up, 'b-*')
plot(B, shift_down, 'r-*')
hold off
xlabel('B /gauss'), ylabel('Abs AC Stark shift /Hz')
legend('Shift on up', 'Shfit on down')
shg