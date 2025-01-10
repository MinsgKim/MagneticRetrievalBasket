classdef External_Magnet2

    properties (Access = private)

        Br = 1.0; % remanence [T], Actually, N35's remanence is 1.17 T
        volume = (0.021)^2*(0.01); % 21x21x10 mm volume
        mu = 4*pi*1e-07; % vacuum permeability, [T*m^2/A]

    end

    methods

        function B_Field = Cal_B_Field(obj, r_vec)

            m = obj.Br *obj.volume / obj.mu;

            r_vec = [r_vec; 0];
            r_norm = norm(r_vec);
            r_hat = r_vec/r_norm;

            B_Field = obj.mu/(4*pi*r_norm^3)*(3*(r_hat*r_hat')-eye(3,3))*[-m; 0; 0];

        end

    end

end