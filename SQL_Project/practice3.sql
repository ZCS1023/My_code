-- 表结构
CREATE TABLE user_login(
    user_id INT,
    login_date DATE
);

-- 用户首次登陆时间表
WITH user_first_login AS(
    SELECT
        user_id,
        MIN(login_date) AS first_login_date
    FROM user_login
    GROUP BY user_id
),

--计算用户登陆时间与首次登陆时间之间的差值
user_login_diff AS(
    SELECT
        a.user_id,
        a.login_date,
        b.first_login_date
        DATEDIFF(a.login_date, b.first_login_date) AS login_diff
    FROM  user_login AS a
    JOIN user_first_login AS b ON a.user_id = b.user_id
)

--计算次日/n日留存率
SELECT
    user_id,
    first_login_date AS 首次登陆时间,
    COUNT(DISTINCT user_id) AS 新增用户数,
    COUNT(DISTINCT CASE WHEN login_diff = 1 THEN user_id ELSE NULL END) AS 留存客户数
    ROUND(
        COUNT(DISTINCT CASE WHEN login_diff = 1 THEN user_id ELSE NULL END) /
        COUNT(DISTINCT user_id) * 100 , 2
    ) AS 次日留存率
FROM user_login_diff
GROUP BY first_login_date
