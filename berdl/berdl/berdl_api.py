import requests


class BERDLAPI:

    def __init__(self, token):
        self.token = token
        self.mcp_url = "https://hub.berdl.kbase.us/apis/mcp/delta/tables/query"
        self.headers = {
            "accept": "application/json",
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
            "User-Agent": "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"
        }
        self.timeout = 1200

    def sql(self, sql, page_limit, offset):
        payload = {
                "query": sql,
                "limit": page_limit,
                "offset": offset
            }
        response = requests.post(
            self.mcp_url,
            headers=self.headers,
            json=payload,
            timeout=self.timeout
        )

        if response.status_code == 200:
            data = response.json()
            results = data.get('result', [])
            return results
        return response

    def get_genome_proteins(self):
        pass

    def get_genome_features(self):
        pass
