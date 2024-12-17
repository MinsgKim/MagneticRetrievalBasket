classdef External_Magnet

    properties

        volume = [];
        remanence = [];
        mu = 4*pi*1e-07;

    end

    methods
        function B = Cal_B(obj, r_vec, m_ext)   % calculate B vector
                % dipole 모델을 사용하여 외부 자기장 계산
                % r_vec: 위치 벡터 (slave magnet 위치 - 외부 자석 위치)
                % m_ext: 외부 자석의 자기 모멘트 (스칼라 값)

            r_norm = norm(r_vec);

            % 외부 자석의 자기 모멘트 방향 (y축 방향으로 가정, 아래쪽을 향함)
            m_vec = [0; -m_ext; 0];

            % 위치 벡터를 3D로 확장
            r_vec_3D = [r_vec; 0];

            % 자기장 계산
            B_full = (obj.mu / (4 * pi)) * ( (3 * r_vec_3D * (dot(m_vec, r_vec_3D)) / r_norm^5) - (m_vec / r_norm^3) );

            % 2D 평면에서의 자기장 (x, y 성분)
            B = B_full(1:2);
        end
    end

end