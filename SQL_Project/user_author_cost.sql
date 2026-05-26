--表结构(四列：打赏日期，用户id，主播id，打赏金额)
CREATE TABLE user_author_live_amt(
    cost_date DATE,  -- 打赏日期
    user_id INT,     -- 用户ID
    author_id INT,   -- 主播ID
    cost_amt INT     -- 打赏金额
)

with t1 as(
    -- 筛选累计打赏超过10万元的用户
    select
        user_id,
        sum(cost_amt) as total_cost  
    from user_author_live_amt  
    group by user_id
    having total_cost > 100000
),

t2 as(
    -- 计算这些用户对每个主播的总打赏金额，并按用户分组排名
    select
        a.user_id,
        a.author_id,
        sum(a.cost_amt) as author_cost,  
        row_number()over(
            partition by a.user_id
            order by sum(a.cost_amt) desc  -- 排序依据需与聚合结果一致
        ) as rank
    from user_author_live_amt as a  -- 正确表名
    join t1 as b on a.user_id = b.user_id
    group by a.user_id,a.author_id
)

-- 输出每个用户的前3主播（包含user_id以区分用户）
select 
    user_id,  -- 保留用户ID，否则无法对应到具体用户
    author_id, 
    author_cost
from t2
where rank <=3
order by user_id, rank;  -- 按用户+排名排序，更清晰