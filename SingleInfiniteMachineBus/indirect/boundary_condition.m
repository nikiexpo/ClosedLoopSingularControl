function res = boundary_condition(x0, xf)

res = [x0(1) - 1.5
        x0(2) - 15
        xf(3) - 0
        xf(4) - 0];
end