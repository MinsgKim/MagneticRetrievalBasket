classdef External_Magnet

    properties

        Br;
        volume;
        mu = 4*pi*1e-07;    % vacuum permeability
        m;

    end

    methods

        function obj = External_Magnet()

            % parameters of an external magnet
            obj.Br = 1.00; % remanence [T], original 1.22 T
            obj.volume = (0.03)^2*0.01; % Volume (3cm x 3cm x 1cm cube)
            obj.m = obj.Br * obj.volume / obj.mu; % magnetic moment

        end


        function B_full = Cal_B(obj, r_vec)

            % parameters of an external magnet
            obj.Br = 1.22; % remanence [T]
            obj.volume = (0.03)^3; % Volume (3cm x 3cm x 3cm cube)
            obj.m = obj.Br * obj.volume / obj.mu; % magnetic moment


            % Dipole 모델을 사용하여 외부 자기장을 계산
            % r_vec: 위치 벡터 (slave magnet 위치 - 외부 자석 위치)
            % m_ext: 외부 자석의 자기 모멘트 (스칼라 값)

            % Calculate the norm of the position vector
            r_norm = norm(r_vec);

            % External magnet's magnetic moment direction (assumed along the y-axis)
            m_vec = [0; -obj.m; 0];

            % Extend r_vec to 3D
            if length(r_vec) == 2
                r_vec_3D = [r_vec; 0];
            else
                r_vec_3D = r_vec;
            end

            % Calculate magnetic field using dipole formula
            B_full = (obj.mu / (4 * pi)) * ( (3 * r_vec_3D * (dot(m_vec, r_vec_3D)) / r_norm^5) - (m_vec / r_norm^3) );

        end



    end

end