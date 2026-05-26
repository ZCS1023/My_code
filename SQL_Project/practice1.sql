-- 2025/08/07 重复练习
-- CREATE TABLE user_author_live_amt(
--     cost_date DATE,
--     user_id INT,
--     author_id INT,
--     cost_amt INT
-- );

WITH t1 AS (
    SELECT
        user_id,
        SUM(cost_amt) AS total_cost
    FROM user_author_live_amt
    GROUP BY user_id
    HAVING total_cost > 100000
),

t2 AS(
    SELECT
        a.user_id,
        a.author_id,
        SUM(a.cost_amt) AS total_author_amt,
        ROW_NUMBER()OVER(
            PARTITION BY a.user_id
            ORDER BY SUM(a.cost_amt) DESC
        ) as author_rank
    FROM user_author_live_amt AS a
    JOIN t1 as b ON a.user_id = b.user_id
    GROUP BY a.user_id, a.author_id
)

SELECT user_id, author_id, total_author_amt
FROM t2
WHERE author_rank <= 3
ORDER BY user_id, author_rank;