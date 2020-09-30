classdef ALLFUNCS
    % main class with static methods with 6 sub-functions
    methods(Static)
        
        function [x_m, fxm] = IntervalHalving(obj_fun, xL, xU, e)
            x_m = (xL + xU)/2; % mean value
            L = abs(xL - xU); % initial value of L
            fxm = obj_fun(x_m); % function value at x_m
            i = 0; % number of interations required to reach an optimum point

            while (L >= e) % looping until we reach the stopping criterion
                L = xU - xL;
                x1 = xL + L/4;
                x2 = xU - L/4; 
                fx1 = obj_fun(x1); % fval at x1
                fx2 = obj_fun(x2); % fval at x2
                i = i + 1; % incrementing the iteration number
                if (fx1 < fxm)
                    xU = x_m;
                    x_m = x1;
                    fxm = fx1;
                elseif (fx2 < fxm)
                    xL = x_m;
                    x_m = x2;
                    fxm = fx2;
                else
                    xL = x1;
                    xU = x2;
                end
            end
        end
        
        function [xe, fxe] = GoldenSection(obj_fun, xL, xU, e)
            i = 0; % counter of iterations
            fL = obj_fun(xL);
            fU = obj_fun(xU);
            R = 0.5*(sqrt(5) - 1); % our golden ratio val
            L = R*(xU - xL);
            x1 = xU - L;
            x2 = xL + L;
            f1 = obj_fun(x1);
            f2 = obj_fun(x2);
            err = inf; % initial value of error span

            while (err >= e)
               i = i + 1; % incrementing counter 
               if (f1 > f2)
                   % changing upper bound
                   xU = x2;
                   fU = f2;
                   % new x2
                   x2 = x1; 
                   f2 = f1;
                   % new x1
                   L = R*(xU - xL);
                   x1 = xU - L;
                   % new f1
                   f1 = obj_fun(x1);
               elseif (f1 < f2)
                   % changing lower bound
                   xL = x1;
                   fL = f1;
                   % new x1
                   x1 = x2;
                   f1 = f2;
                   % new x2
                   L = R*(xU - xL);
                   x2 = xL + L;
                   % new f2
                   f2 = obj_fun(x2);
               else
                   xU = (x1 + x2)/2;
                   xL = xU; % to break out of loop
               end
               err = abs(xU - xL); % updating the span
            end

            xe = (x1 + x2)/2; % the estimated optimum point
            fxe = obj_fun(xe); % the minimum fval estimated
        end
        
        function [xm, fxm] = QuadraticEstimation(obj_fun, xL, xU)
            x1 = xL; % initial val of x1
            x2 = (xL + xU)/2; % initial val of x2
            x3 = xU; % initial val of x3
            
            fx1 = obj_fun(x1);
            fx2 = obj_fun(x2);
            fx3 = obj_fun(x3);

            % determining constants of quadratic equation
            a0 = fx1; 
            a1 = (fx2 - fx1)/(x2 - x1);
            a2 = 1/(x3 - x2) * ((fx3 - fx1)/(x3 - x1) - (fx2 - fx1)/(x2 - x1));
            xm = (x2 + x1)/2 - a1/(2*a2); % estimated optimum point
            fxm = obj_fun(xm);
        end
        
        function [x0, fx0] = NewtonRaphson(obj_fun, x, x0, e)
            i = 0; % iteration counter
            f_der = diff(obj_fun); % derivative of obj_fun
            err = 10000; % initial value for difference error
            
            for j = 1:100 % looping, might be increased if necessary
               i = i + 1; % incrementing counter
               fx0 = vpa(subs(obj_fun,x,x0)); % val of obj_fun at x0 by replacing x
               f_derx0 = vpa(subs(f_der,x,x0)); % val of f_der at x0 by replacing x
               new_x = x0 - fx0/f_derx0; % new point
               err = abs(new_x - x0); % difference error
               if err < e
                   break % breaking the loop
               end
               x0 = new_x;
            end
        end
        
        function [z, fz] = Bisection(obj_fun, x, xL, xU, e)
            i = 0; % iteration counter
            f_der = diff(obj_fun); % derivative of obj_fun

            for j = 1:100 % looping, might be increased if necessary
                i = i + 1; % incrementing counter
                z = (xL + xU)/2; % bisection value
                f_derz = vpa(subs(f_der,x,z)); % value of d_der at z via substituting x
                if abs(f_derz) <= e % checking for stopping criterion
                   break
                elseif f_derz < 0 % min on the right side
                    xL = z;
                elseif f_derz > 0 % min on the left side
                    xU = z;
                end
            end
            fz = vpa(subs(obj_fun,x,z)); % value of obj_fun at z via substituting x
        end
        
        function [z, fz] = Secant(obj_fun, x, xL, xU, e)
            i = 0; % iteration counter
            f_der = diff(obj_fun); % derivative of obj_fun

            for j = 1:100 % looping, might be increased if necessary
               i = i + 1; % incrementing counter
               f_derxL = vpa(subs(f_der,x,xL)); % value of d_der at xL via substituting x
               f_derxU = vpa(subs(f_der,x,xU)); % value of d_der at xU via substituting x
               z = xU - f_derxU/(f_derxU - f_derxL/(xU - xL)); % secant value
               f_derz = vpa(subs(f_der,x,z)); % value of d_der at z via substituting x
               if abs(f_derz) <= e % checking for stopping criterion
                   break
               elseif f_derz < 0 % min on the right side
                   xL = z;
               elseif f_derz > 0 % min on the left side
                   xU = z;
               end
            end
            fz = vpa(subs(obj_fun,x,z)); % value of obj_fun at z via substituting x
        end
        
    end
end