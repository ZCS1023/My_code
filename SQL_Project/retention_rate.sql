--计算首次登录日期
WITH user_first_login AS(
    SELECT
        user_id,
        MIN(login_date) AS first_login
    FROM user_login
    GROUP BY user_id
),

--计算当前登录日期和首次登录日期之间的差
user_login_with_diff AS(
    SELECT
        a.user_id,
        a.login_date,
        b.first_login,
        DATEDIFF(a.login_date, b.first_login) AS day_diff
    FROM user_login a
    JOIN user_first_login ON a.user_id = b.user_id
)

--以首次登陆日期字段来分组，统计day_diff = 1（次日）的留存率情况
SELECT
    first_login AS 首次登录时间,
    COUNT(DISTINCT user_id) AS 新增用户数,
    COUNT(DISTINCT CASE WHEN day_diff = 1 THEN user_id ELSE NULL END) AS 次日留存用户数,
    ROUND(
        COUNT(DISTINCT CASE WHEN day_diff = 1 THEN user_id ELSE NULL END) 
        / COUNT(DISTINCT user_id) *100, 2
    ) AS 次日留存率
FROM user_login_with_diff
GROUP BY first_login
ORDER BY first_login ASC