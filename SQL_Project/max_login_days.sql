-- 题目：求用户最大连续登陆天数
-- CREATE TABLE user_login (
--     user_id INT,       -- 用户ID
--     login_date DATE    -- 登录日期
-- );

-- 查询每个用户的连续最大登录天数

WITH t1 AS(
    SELECT
        DISTINCT user_id, login_date
    FROM user_login
),

t2 AS(
    SELECT
        user_id,
        login_date,
        row_number()OVER(
            PARTITION BY user_id,
            ORDER BY login_date
        ) as date_rank
    FROM t1
),

t3 AS(
    SELECT
        user_id,
        login_date,
        DATE_SUB(login_date, INTERVAL date_rank DAY) AS group_date
    FROM t2
),

t4 AS(
    SELECT
        user_id,
        group_date,
        COUNT(*) AS continuous_days
    FROM t3
    GROUP BY user_id, group_date
)

SELECT
    user_id,
    MAX(continuous_days) AS max_login_days
FROM t4
GROUP BY user_id